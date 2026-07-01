use crate::alignment;
use crate::g2t::{G2TTree, GuideExon};
use crate::sw::{self, Anchor};
use crate::types::{ReadId, RefId, Tid};
use anyhow::Result;
use ahash::RandomState;
use crate::types::{HashMap, HashMapExt, HashSet, HashSetExt};
use noodles::sam::alignment::record::cigar::op::Kind as CigarKind;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(clippy::enum_variant_names)]
pub enum ExonStatus {
    FirstExon = 0,
    MiddleExon = 1,
    LastExon = 2,
    OnlyExon = 3,
    InsExon = 4,
    GapExon = 5,
    LeftClipExon = 6,
    RightClipExon = 7,
}

#[derive(Debug, Clone)]
pub struct ReadEvaluationConfig {
    pub max_clip: u32,
    pub max_ins: u32,
    pub max_junc_gap: u32,
    pub similarity_threshold: f32,
    pub ignore_small_exons: bool,
    pub small_exon_size: u32,
    #[allow(dead_code)]
    pub print: bool,
    #[allow(dead_code)]
    pub name: Vec<u8>,
    pub soft_clips: bool,
    pub strict: bool,
    pub use_fasta: bool,
    /// Whether to apply the similarity filter (retain matches by
    /// `total_coverage/total_operations > threshold` and set
    /// `similarity_score = x²·(junc_hits+1)`). Mirrors C++
    /// `config.filter_by_similarity`, which is `(similarity_threshold < 1.0)`:
    /// short reads use threshold 1.0 → filter DISABLED (keep all matches, leave
    /// `similarity_score` at its default), so short-read projection is not
    /// re-weighted by junc_hits (C++ only rewrites AS for long reads).
    pub filter_by_similarity: bool,
    /// Multiplicative discount in `(0, 1]` applied per internal junction mismatch
    /// (`junc_misses`): `similarity_score *= junc_miss_discount^junc_misses`.
    /// `1.0` = off (default, original behavior); smaller values sharpen isoform
    /// discrimination by penalizing transcripts the read's junctions disagree with.
    pub junc_miss_discount: f64,
}

impl ReadEvaluationConfig {
    /// Parameters matching C++ ShortReadEvaluator.
    pub fn short_read() -> Self {
        Self {
            max_clip: 5,
            max_ins: 0,
            max_junc_gap: 0,
            // C++ short-read `similarity_threshold` = 1.0, which sets
            // `filter_by_similarity = (threshold < 1.0) = false`: the similarity
            // filter is DISABLED for short reads (all matches kept, no junc_hits
            // re-weighting). This is a SENTINEL, not a "require perfect match"
            // gate — the filter is skipped entirely rather than dropping reads.
            similarity_threshold: 1.0,
            filter_by_similarity: false,
            ignore_small_exons: false,
            small_exon_size: 0,
            print: false,
            name: Vec::new(),
            soft_clips: true,
            strict: false,
            use_fasta: false,
            junc_miss_discount: 1.0,
        }
    }

    /// Parameters matching C++ LongReadEvaluator.
    #[allow(dead_code)]
    pub fn long_read() -> Self {
        Self {
            max_clip: 40,
            max_ins: 40,
            max_junc_gap: 40,
            similarity_threshold: 0.60,
            filter_by_similarity: true,
            ignore_small_exons: true,
            small_exon_size: 35,
            print: false,
            name: Vec::new(),
            soft_clips: true,
            strict: false,
            use_fasta: false,
            junc_miss_discount: 1.0,
        }
    }
}

impl Default for ReadEvaluationConfig {
    fn default() -> Self {
        Self::short_read()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum CigarOp {
    #[default]
    Match,
    Ins,
    Del,
    #[allow(dead_code)]
    RefSkip,
    SoftClip,
    #[allow(dead_code)]
    HardClip,
    #[allow(dead_code)]
    Pad,
    #[allow(dead_code)]
    Equal,
    #[allow(dead_code)]
    Diff,
    MatchOverride,
    DelOverride,
    InsOverride,
    ClipOverride,
}

#[derive(Debug, Clone, Default)]
pub struct Cigar {
    pub ops: Vec<(u32, CigarOp)>,
}

impl Cigar {
    pub fn add_operation(&mut self, len: u32, op: CigarOp) {
        if len == 0 {
            return;
        }
        if let Some((prev_len, prev_op)) = self.ops.last_mut()
            && *prev_op == op
        {
            *prev_len += len;
            return;
        }
        self.ops.push((len, op));
    }

    #[allow(dead_code)]
    pub fn prepend_operation(&mut self, len: u32, op: CigarOp) {
        if len == 0 {
            return;
        }
        if let Some((prev_len, prev_op)) = self.ops.first_mut()
            && *prev_op == op
        {
            *prev_len += len;
            return;
        }
        self.ops.insert(0, (len, op));
    }

    #[allow(dead_code)]
    pub fn reverse(&mut self) {
        self.ops.reverse();
    }

    /// Number of query (read) bases consumed by this CIGAR, i.e. the sum of
    /// lengths of query-consuming operations (M/I/=/X and their rescue
    /// overrides). Soft/hard clips and padding are excluded, so this is the
    /// length of the aligned portion of the read.
    pub fn query_consumed(&self) -> u32 {
        self.ops
            .iter()
            .map(|(len, op)| match op {
                CigarOp::Match
                | CigarOp::Equal
                | CigarOp::Diff
                | CigarOp::MatchOverride
                | CigarOp::Ins
                | CigarOp::InsOverride => *len,
                _ => 0,
            })
            .sum()
    }

    #[allow(dead_code)]
    #[allow(clippy::inherent_to_string)]
    pub fn to_string(&self) -> String {
        if self.ops.is_empty() {
            return "*".to_string();
        }
        let mut out = String::new();
        for (len, op) in &self.ops {
            out.push_str(&len.to_string());
            out.push(match op {
                CigarOp::Match => 'M',
                CigarOp::Ins => 'I',
                CigarOp::Del => 'D',
                CigarOp::RefSkip => 'N',
                CigarOp::SoftClip => 'S',
                CigarOp::HardClip => 'H',
                CigarOp::Pad => 'P',
                CigarOp::Equal => '=',
                CigarOp::Diff => 'X',
                CigarOp::MatchOverride => ',',
                CigarOp::DelOverride => '.',
                CigarOp::InsOverride => '/',
                CigarOp::ClipOverride => ';',
            });
        }
        out
    }
}

#[derive(Debug, Clone, Default)]
pub struct AlignInfo {
    pub fwpos: u32,
    pub rcpos: u32,
    pub strand: char,
    pub cigar: Option<Cigar>,
    pub is_reverse: bool,
    pub primary_alignment: bool,
    pub clip_score: i32,
    pub similarity_score: f64,
    pub hit_index: i32,
}

#[derive(Debug, Clone, Default)]
pub struct ExonChainMatch {
    #[allow(dead_code)]
    pub created: bool,
    pub align: AlignInfo,
    pub total_coverage: f64,
    pub total_operations: f64,
    pub ref_consumed: i32,
    pub junc_hits: i32,
    /// Number of internal exon boundaries where the read does NOT match this
    /// transcript exactly (a tolerated gap/insertion at an internal junction) —
    /// evidence the read came from a different isoform. Used to discount the
    /// similarity score (see `junc_miss_discount`).
    pub junc_misses: i32,
    pub transcript_len: i32,
    pub prev_op: CigarOp,
}

#[derive(Debug, Clone)]
pub struct EvalSegment {
    #[allow(dead_code)]
    pub is_valid: bool,
    pub has_gexon: bool,
    #[allow(dead_code)]
    pub has_qexon: bool,
    pub gexon: Option<GuideExon>,
    pub qexon: alignment::Segment,
    pub status: ExonStatus,
    #[allow(dead_code)]
    pub is_small_exon: bool,
    pub cigar: Cigar,
    pub score: i32,
}

impl Default for EvalSegment {
    fn default() -> Self {
        Self {
            is_valid: false,
            has_gexon: false,
            has_qexon: false,
            gexon: None,
            qexon: alignment::Segment { start: 0, end: 0 },
            status: ExonStatus::OnlyExon,
            is_small_exon: false,
            cigar: Cigar::default(),
            score: 0,
        }
    }
}

#[derive(Debug, Default)]
pub struct TidData {
    pub elim: bool,
    pub segments: Vec<EvalSegment>,
    pub match_info: ExonChainMatch,
    pub has_left_clip: bool,
    pub has_right_clip: bool,
}

#[derive(Debug)]
pub struct ReadAln {
    pub strand: char,
    pub strand_from_tag: bool,
    pub refid: RefId,
    #[allow(dead_code)]
    pub nh: u16,
    pub segs: Vec<alignment::Segment>,
    /// Raw CIGAR ops: (length, noodles CigarKind).
    pub cigar_ops: Vec<(u32, CigarKind)>,
    /// Query sequence bytes (optional; needed for clip rescue).
    pub sequence: Option<Vec<u8>>,
    /// Read name (used for deterministic tie-breaking in multi-hit selection).
    pub name: Vec<u8>,
}

/// Reusable scratch space for one evaluation pass.
/// Allocate once per worker thread; cleared between reads automatically.
#[doc(hidden)]
pub struct EvalContext {
    /// Final output: transcript → chain match.  Filled by evaluate(); drained by caller.
    pub matches: HashMap<Tid, ExonChainMatch>,
    /// Per-strand working map: tid → accumulated TidData.
    pub(crate) data: HashMap<Tid, TidData>,
    /// Working set of matched tids for the current exon.
    pub(crate) candidate_tids: HashSet<Tid>,
    /// Scratch buffer for guide exon queries.
    pub(crate) guide_exons: HashMap<Tid, GuideExon>,
    /// Reusable scratch buffer for tid keys during evaluate_exon_chains.
    tids_buf: Vec<Tid>,
    /// Reusable ksw2 aligner (workspace + result buffers) for clip rescue SW.
    pub(crate) aligner: ksw2rs::Aligner,
    /// Reusable scratch buffers for SW encoding.
    pub(crate) sw_bufs: sw::SwBufs,
}

impl EvalContext {
    pub fn new() -> Self {
        Self {
            matches: HashMap::new(),
            data: HashMap::new(),
            candidate_tids: HashSet::new(),
            guide_exons: HashMap::new(),
            tids_buf: Vec::new(),
            aligner: ksw2rs::Aligner::new(),
            sw_bufs: sw::SwBufs::new(),
        }
    }
}

pub fn get_exon_status(exon_count: usize, j: usize) -> ExonStatus {
    if exon_count == 1 {
        ExonStatus::OnlyExon
    } else if j == 0 {
        ExonStatus::FirstExon
    } else if j < (exon_count - 1) {
        ExonStatus::MiddleExon
    } else {
        ExonStatus::LastExon
    }
}

pub fn get_strands_to_check(read: &ReadAln, long_reads: bool) -> Vec<char> {
    // C++ always checks both strands for long reads, regardless of splice-strand
    // tags. The strand tag is an informed guess that may be wrong for long reads,
    // so both are always evaluated and filtered by similarity.
    if long_reads {
        return vec!['+', '-'];
    }
    if read.strand == '+' {
        return vec!['+'];
    }
    if read.strand == '-' {
        return vec!['-'];
    }
    vec!['+', '-']
}

pub fn get_clips(
    read: &ReadAln,
    config: &ReadEvaluationConfig,
) -> Result<(bool, bool, u32, u32)> {
    let mut has_left_clip = false;
    let mut has_right_clip = false;
    let mut n_left_clip = 0;
    let mut n_right_clip = 0;

    let ops = &read.cigar_ops;

    if ops.is_empty() {
        return Ok((false, false, 0, 0));
    }

    let (left0_len, left0_kind) = ops.first().unwrap();
    let left1 = ops.get(1);
    // C++ sets has_left_clip = USE_FASTA (only true when FASTA mode is on).
    // Without FASTA, soft clips are detected (n_left_clip is set) but
    // has_left_clip stays false, which changes how build_cigar_match handles
    // boundary insertions (SoftClip vs Ins).
    if *left0_kind == CigarKind::HardClip {
        if let Some((op_len, op_kind)) = left1
            && *op_kind == CigarKind::SoftClip
        {
            has_left_clip = config.use_fasta;
            n_left_clip = *op_len;
        }
    } else if *left0_kind == CigarKind::SoftClip {
        has_left_clip = config.use_fasta;
        n_left_clip = *left0_len;
    }

    let (right0_len, right0_kind) = ops.last().unwrap();
    let right1 = if ops.len() >= 2 { ops.get(ops.len() - 2) } else { None };
    if *right0_kind == CigarKind::HardClip {
        if let Some((op_len, op_kind)) = right1
            && *op_kind == CigarKind::SoftClip
        {
            has_right_clip = config.use_fasta;
            n_right_clip = *op_len;
        }
    } else if *right0_kind == CigarKind::SoftClip {
        has_right_clip = config.use_fasta;
        n_right_clip = *right0_len;
    }

    Ok((has_left_clip, has_right_clip, n_left_clip, n_right_clip))
}

#[allow(clippy::too_many_arguments)]
pub fn get_intervals(
    ctx: &mut EvalContext,
    read: &ReadAln,
    j: usize,
    exon_count: usize,
    config: &ReadEvaluationConfig,
    g2t: &G2TTree,
    refid: RefId,
    strand: char,
    has_left_clip: bool,
    has_right_clip: bool,
    start_ignore: &mut EvalSegment,
    end_ignore: &mut EvalSegment,
    has_gap: &mut bool,
    failure: &mut bool,
) {
    let _ = (start_ignore, end_ignore);
    let qexon = read.segs[j];
    let status = get_exon_status(exon_count, j);
    let is_small_exon = (qexon.end - qexon.start) <= config.small_exon_size;
    let data_empty = ctx.data.is_empty();

    g2t.get_guide_exons(refid, strand, qexon.start, qexon.end, config, status, &mut ctx.guide_exons);
    if !ctx.guide_exons.is_empty() {
        ctx.candidate_tids.clear();

        for (tid, gexon) in ctx.guide_exons.drain() {
            ctx.candidate_tids.insert(tid);

            if data_empty && status != ExonStatus::InsExon {
                let tid_data = TidData {
                    has_left_clip,
                    has_right_clip,
                    ..Default::default()
                };
                if gexon.left_gap > 0 {
                    *has_gap = true;
                }
                ctx.data.insert(tid, tid_data);
            } else if !ctx.data.contains_key(&tid) || ctx.data.get(&tid).map(|d| d.elim).unwrap_or(true) {
                continue;
            }

            if matches!(status, ExonStatus::LastExon | ExonStatus::OnlyExon)
                && gexon.right_gap > 0
            {
                *has_gap = true;
            }

            if let Some(tid_data) = ctx.data.get_mut(&tid) {
                let seg = EvalSegment {
                    is_valid: true,
                    has_qexon: true,
                    has_gexon: true,
                    gexon: Some(gexon),
                    qexon,
                    status,
                    is_small_exon,
                    ..Default::default()
                };
                tid_data.segments.push(seg);
            }
        }

        for (tid, td) in ctx.data.iter_mut() {
            if !ctx.candidate_tids.contains(tid) {
                td.elim = true;
            }
        }

        return;
    }

    // guide exons empty
    if status != ExonStatus::OnlyExon && config.ignore_small_exons && is_small_exon {
        if status == ExonStatus::MiddleExon {
            if ctx.data.is_empty() {
                *failure = true;
                return;
            }
            for td in ctx.data.values_mut() {
                let seg = EvalSegment {
                    is_valid: true,
                    has_qexon: true,
                    has_gexon: false,
                    gexon: None,
                    qexon,
                    status: ExonStatus::InsExon,
                    is_small_exon: true,
                    ..Default::default()
                };
                td.segments.push(seg);
            }
        } else {
            *failure = true;
        }
        return;
    }

    *failure = true;
}

fn seq_slice_from_shared(shared: Option<&[u8]>, read: &ReadAln, start: usize, len: usize) -> Vec<u8> {
    if let Some(seq) = shared
        && !seq.is_empty()
    {
        let end = (start + len).min(seq.len());
        return seq[start..end].to_vec();
    }

    if let Some(seq) = &read.sequence {
        let end = (start + len).min(seq.len());
        return seq[start..end].to_vec();
    }
    Vec::new()
}

#[allow(clippy::too_many_arguments)]
pub fn left_clip_rescue(
    tid_data: &mut TidData,
    strand: char,
    g2t: &G2TTree,
    refid: RefId,
    tid: Tid,
    n_left_clip: u32,
    config: &ReadEvaluationConfig,
    read: &ReadAln,
    shared_seq: Option<&[u8]>,
    aligner: &mut ksw2rs::Aligner,
    sw_bufs: &mut sw::SwBufs,
) {
    if !config.use_fasta {
        return;
    }
    if tid_data.segments.is_empty() {
        return;
    }

    let segment = &mut tid_data.segments[0];
    if !segment.has_gexon {
        tid_data.has_left_clip = false;
        return;
    }
    let gexon = segment.gexon.clone().unwrap();
    if gexon.left_gap > 0 {
        tid_data.has_left_clip = false;
        return;
    }

    let total_bases = (n_left_clip as i32 + gexon.left_ins).max(0) as usize;
    if total_bases == 0 {
        tid_data.has_left_clip = false;
        return;
    }

    let mut remaining_qseq = seq_slice_from_shared(shared_seq, read, 0, total_bases);
    if remaining_qseq.is_empty() {
        tid_data.has_left_clip = false;
        return;
    }

    let current_gexon = gexon.clone();
    // Multi-exon clip rescue: collect neighboring exons until guide >= query length.
    // If multi-exon alignment fails, fall back to single-exon (baseline behavior).
    {
        if !remaining_qseq.is_empty() {
            let mut gseq: Vec<u8> = Vec::new();
            let mut curr = current_gexon.clone();
            let mut first_gexon: Option<GuideExon> = None;
            let mut n_collected = 0u32;

            while remaining_qseq.len() > gseq.len() {
                let has_neighbor = if strand == '+' { curr.has_prev } else { curr.has_next };
                if !has_neighbor {
                    if n_collected == 0 { tid_data.has_left_clip = false; return; }
                    break;
                }

                let next = if strand == '+' {
                    g2t.get_guide_exon_for_tid(refid, strand, tid, curr.prev_start, curr.prev_end)
                } else {
                    g2t.get_guide_exon_for_tid(refid, strand, tid, curr.next_start, curr.next_end)
                };
                let Some(next) = next else {
                    if n_collected == 0 { tid_data.has_left_clip = false; return; }
                    break;
                };
                if next.seq.is_empty() {
                    if n_collected == 0 { tid_data.has_left_clip = false; return; }
                    break;
                }

                // Prepend: new exon sequence goes before accumulated guide
                let mut new_gseq = Vec::with_capacity(next.seq.len() + gseq.len());
                new_gseq.extend_from_slice(&next.seq);
                new_gseq.extend_from_slice(&gseq);
                gseq = new_gseq;
                n_collected += 1;
                if first_gexon.is_none() {
                    first_gexon = Some(next.clone());
                }
                curr = next;
            }

            if gseq.is_empty() {
                tid_data.has_left_clip = false;
                return;
            }

            // Trim guide: keep the rightmost query_len + 40 bytes (matching C++ align_reversed windowing).
            let window = remaining_qseq.len() + 40;
            if gseq.len() > window {
                let start = gseq.len() - window;
                gseq.drain(..start);
            }

            let result = sw::smith_waterman(aligner, sw_bufs, &remaining_qseq, &gseq, Anchor::End);
            if result.score < 10 || result.zdropped {
                tid_data.has_left_clip = false;
                return;
            }

            let mut query_consumed: i32 = 0;
            let mut ref_consumed: i32 = 0;
            for (len, op) in &result.cigar.ops {
                let len = *len as i32;
                match op {
                    CigarOp::MatchOverride | CigarOp::InsOverride | CigarOp::ClipOverride => {
                        query_consumed += len;
                    }
                    _ => {}
                }
                match op {
                    CigarOp::MatchOverride | CigarOp::DelOverride => {
                        ref_consumed += len;
                    }
                    _ => {}
                }
            }

            let q_len = remaining_qseq.len() as i32;
            let _left_clip_bases = (q_len - query_consumed).max(0) as u32;

            let mut dummy_gexon = gexon.clone();
            dummy_gexon.start = gexon.start.saturating_sub(ref_consumed as u32);
            dummy_gexon.end = gexon.start;
            dummy_gexon.pos = gexon.pos_start.saturating_sub(ref_consumed as u32);

            // Build the CIGAR: prepend soft-clip for unaligned bases, then
            // iterate the reversed ksw2 CIGAR (C++ iterates n_cigar-1 down to 0).
            let mut cigar = Cigar::default();

            // in C++
            // if left_clip_bases > 0 {
            //     cigar.ops.push((left_clip_bases, CigarOp::ClipOverride))
            // }

            if let Some(&(len, ref op)) = result.cigar.ops.first() {
                match op {
                    CigarOp::DelOverride => { /* discard leading D*/ }
                    CigarOp::InsOverride => {
                        cigar.ops.push((len, CigarOp::ClipOverride));
                    }
                    _ => {
                        cigar.ops.push((len, op.clone()))
                    }
                }
                for (len, op) in result.cigar.ops.iter().skip(1) {
                    cigar.ops.push((*len, op.clone()))
                }
            } // else: empty CIGAR, nothing to append

            let left_clip = EvalSegment {
                is_valid: true,
                has_qexon: false,
                has_gexon: true,
                gexon: Some(dummy_gexon),
                qexon: alignment::Segment::default(), // not used
                status: ExonStatus::LeftClipExon,
                is_small_exon: remaining_qseq.len() as u32 <= config.small_exon_size,
                cigar,
                score: result.score,
            };

            tid_data.segments.insert(0, left_clip);
            // Match C++: zero out left_ins on the original first segment (now at index 1)
            // so build_cigar_match sees left_ins=0 and increments junc_hits.
            if let Some(seg) = tid_data.segments.get_mut(1)
                && let Some(ref mut ge) = seg.gexon
                && ge.left_ins > 0
            {
                ge.left_ins = 0;
            }
            tid_data.has_left_clip = true;

            if result.start_i > 1 {
                remaining_qseq.truncate(result.start_i as usize);
            } else {
                remaining_qseq.clear();
            }
        }
    }

    // Do NOT add ClipOverride for unrescued bases. The unrescued query bases
    // from the leading soft clip will naturally pass through the merge as S
    // ops when they hit the non-override M ops of the ideal (via pad_ideal_for_leading_clips).
    // Adding ClipOverride here would cause wrong ordering (before the rescue match).
}

#[allow(clippy::too_many_arguments)]
pub fn right_clip_rescue(
    tid_data: &mut TidData,
    strand: char,
    g2t: &G2TTree,
    refid: RefId,
    tid: Tid,
    n_right_clip: u32,
    config: &ReadEvaluationConfig,
    read: &ReadAln,
    shared_seq: Option<&[u8]>,
    aligner: &mut ksw2rs::Aligner,
    sw_bufs: &mut sw::SwBufs,
) {
    if !config.use_fasta {
        return;
    }
    if tid_data.segments.is_empty() {
        return;
    }

    let last_idx = tid_data.segments.len() - 1;
    let segment = &mut tid_data.segments[last_idx];
    if !segment.has_gexon {
        tid_data.has_right_clip = false;
        return;
    }
    let gexon = segment.gexon.clone().unwrap();
    if gexon.right_gap > 0 {
        tid_data.has_right_clip = false;
        return;
    }

    let total_bases = (n_right_clip as i32 + gexon.right_ins).max(0) as usize;
    if total_bases == 0 {
        tid_data.has_right_clip = false;
        return;
    }
    let seq_len = if let Some(seq) = shared_seq
        && !seq.is_empty()
    {
        seq.len()
    } else {
        read.sequence.as_ref().map_or(0, |s| s.len())
    };
    let start_pos = seq_len.saturating_sub(total_bases);
    let mut remaining_qseq = seq_slice_from_shared(shared_seq, read, start_pos, total_bases);
    if remaining_qseq.is_empty() {
        tid_data.has_right_clip = false;
        return;
    }

    let current_gexon = gexon.clone();
    // Multi-exon clip rescue: collect neighboring exons until guide >= query length.
    // If multi-exon alignment fails, fall back to single-exon (baseline behavior).
    {
        if !remaining_qseq.is_empty() {
            let mut gseq: Vec<u8> = Vec::new();
            let mut curr = current_gexon.clone();
            let mut first_gexon: Option<GuideExon> = None;
            let mut n_collected = 0u32;

            while remaining_qseq.len() > gseq.len() {
                let has_neighbor = if strand == '+' { curr.has_next } else { curr.has_prev };
                if !has_neighbor {
                    if n_collected == 0 { tid_data.has_right_clip = false; return; }
                    break;
                }

                let next = if strand == '+' {
                    g2t.get_guide_exon_for_tid(refid, strand, tid, curr.next_start, curr.next_end)
                } else {
                    g2t.get_guide_exon_for_tid(refid, strand, tid, curr.prev_start, curr.prev_end)
                };
                let Some(next) = next else {
                    if n_collected == 0 { tid_data.has_right_clip = false; return; }
                    break;
                };
                if next.seq.is_empty() {
                    if n_collected == 0 { tid_data.has_right_clip = false; return; }
                    break;
                }

                gseq.extend_from_slice(&next.seq);
                n_collected += 1;
                if first_gexon.is_none() {
                    first_gexon = Some(next.clone());
                }
                curr = next;
            }

            if gseq.is_empty() {
                tid_data.has_right_clip = false;
                return;
            }

            // Trim guide to query_len + 40 (matching C++).
            let window = remaining_qseq.len() + 40;
            if gseq.len() > window { gseq.truncate(window); }

            let result = sw::smith_waterman(aligner, sw_bufs, &remaining_qseq, &gseq, Anchor::Start);
            if result.score < 10 || result.zdropped {
                tid_data.has_right_clip = false;
                return;
            }

            let mut query_consumed: i32 = 0;
            let mut ref_consumed: i32 = 0;
            for (len, op) in &result.cigar.ops {
                let len = *len as i32;
                match op {
                    CigarOp::MatchOverride | CigarOp::InsOverride | CigarOp::ClipOverride => {
                        query_consumed += len;
                    }
                    _ => {}
                }
                match op {
                    CigarOp::MatchOverride | CigarOp::DelOverride => {
                        ref_consumed += len;
                    }
                    _ => {}
                }
            }

            let q_len = remaining_qseq.len() as i32;
            let _right_clip_bases = (q_len - query_consumed).max(0) as u32;

            let mut dummy_gexon = gexon.clone();
            dummy_gexon.start = gexon.end;
            dummy_gexon.end = gexon.end + ref_consumed as u32;
            dummy_gexon.pos = gexon.pos_start.saturating_sub(ref_consumed as u32);

            let mut cigar = Cigar::default();

            let n = result.cigar.ops.len();
            for (i, (len, op)) in result.cigar.ops.iter().enumerate() {
                if i == n - 1 {
                    match op {
                        CigarOp::DelOverride => { /* trailing D - discard */}
                        CigarOp::InsOverride => {
                            cigar.ops.push((*len, CigarOp::ClipOverride));
                        }
                        _ => {
                            cigar.ops.push((*len, op.clone()));
                        }
                    }
                } else {
                    cigar.ops.push((*len, op.clone()));
                }
            }

            // in C++
            // if right_clip_bases > 0 {
            //     cigar.ops.push((right_clip_bases, CigarOp::ClipOverride))
            // }

            let right_clip = EvalSegment {
                is_valid: true,
                has_qexon: false,
                has_gexon: true,
                gexon: Some(dummy_gexon),
                qexon: alignment::Segment::default(),
                status: ExonStatus::RightClipExon,
                is_small_exon: remaining_qseq.len() as u32 <= config.small_exon_size,
                cigar,
                score: result.score,
            };
            tid_data.segments.push(right_clip);
            // Match C++: zero out right_ins on the original last segment
            // so build_cigar_match sees right_ins=0 and increments junc_hits.
            let orig_last = tid_data.segments.len() - 2;
            if let Some(ref mut ge) = tid_data.segments[orig_last].gexon
                && ge.right_ins > 0
            {
                ge.right_ins = 0;
            }
            tid_data.has_right_clip = true;

            let query_consumed = result.end_i + 1;
            if query_consumed < remaining_qseq.len() as i32 {
                remaining_qseq.drain(..query_consumed as usize);
            } else {
                remaining_qseq.clear();
            }
        }
    }

    // Do NOT add ClipOverride for unrescued bases. The unrescued query bases
    // from the trailing soft clip will naturally pass through the merge as S
    // ops. Adding ClipOverride here causes phantom S ops when real is exhausted
    // but the ideal still has ClipOverride entries.

    // this shouldn't be true; the query only takes from existing soft clips
}

pub fn correct_for_gaps(
    tid_data: &mut TidData,
    tid: Tid,
    config: &ReadEvaluationConfig,
    g2t: &G2TTree,
    strand: char,
    refid: RefId,
    long_reads: bool,
) {
    // Take ownership of segments to avoid cloning each one.
    let segments = std::mem::take(&mut tid_data.segments);
    let mut out: Vec<EvalSegment> = Vec::with_capacity(segments.len() + 4);
    let mut prev_gexon: Option<GuideExon> = None;

    for segment in segments {
        if tid_data.elim {
            return;
        }
        if !segment.has_gexon {
            out.push(segment);
            continue;
        }

        let gexon = segment.gexon.as_ref().unwrap();

        if let Some(ref prev) = prev_gexon {
            let gap = gexon.exon_id.saturating_sub(prev.exon_id);

            // here, to match C++, should only turn on for --lr (not --lr_hq)
            // specifically, when small_exon_size > 0
            if !long_reads {
                if gap != 1 {
                    tid_data.elim = true;
                    return;
                }
            } else {
                if gap > 2 {
                    tid_data.elim = true;
                    return;
                }
                if gap == 2 {
                    let (gap_start, gap_end) = if strand == '+' {
                        (gexon.prev_start, gexon.prev_end)
                    } else {
                        (gexon.next_start, gexon.next_end)
                    };

                    if (gap_start == 0 && gap_end == 0)
                        || gap_end.saturating_sub(gap_start) > config.small_exon_size
                    {
                        tid_data.elim = true;
                        return;
                    }

                    let prev_exon =
                        g2t.get_guide_exon_for_tid(refid, strand, tid, gap_start, gap_end);
                    let Some(prev_exon) = prev_exon else {
                        tid_data.elim = true;
                        return;
                    };

                    let gap_seg = EvalSegment {
                        is_valid: true,
                        has_gexon: true,
                        has_qexon: false,
                        gexon: Some(prev_exon),
                        status: ExonStatus::GapExon,
                        is_small_exon: (gap_end.saturating_sub(gap_start)
                            <= config.small_exon_size),
                        ..Default::default()
                    };

                    out.push(gap_seg);
                }
            }
        }

        prev_gexon = segment.gexon.clone();
        out.push(segment);
    }

    tid_data.segments = out;
}

pub fn create_match(
    gexon: &GuideExon,
    strand: char,
    read_strand: char,
    ins_exon_len: u32,
) -> ExonChainMatch {
    let mut match_info = ExonChainMatch::default();
    match_info.align.fwpos = gexon.pos;
    match_info.align.rcpos = gexon.pos;
    match_info.transcript_len = gexon.transcript_len as i32;
    match_info.align.strand = strand;
    match_info.align.cigar = Some(Cigar::default());
    match_info.align.is_reverse = strand != read_strand;
    match_info.align.similarity_score = 0.0;
    match_info.total_coverage = 0.0;
    match_info.total_operations = ins_exon_len as f64;
    match_info.ref_consumed = 0;
    match_info.prev_op = CigarOp::Match;
    match_info.junc_hits = 0;
    match_info.junc_misses = 0;
    match_info
}

pub fn build_cigar_match(
    segment: &EvalSegment,
    has_left_clip: bool,
    has_right_clip: bool,
    match_info: &mut ExonChainMatch,
    first_match: bool,
    last_match: bool,
) {
    let qstart = segment.qexon.start;
    let qend = segment.qexon.end;
    let gexon = segment.gexon.as_ref().unwrap();
    let gstart = gexon.start;
    let gend = gexon.end;

    // Use left_ins/left_gap/right_ins/right_gap fields from gexon, matching C++.
    // C++ uses these fields (not positional comparison) to decide whether to
    // add ins/del/softclip at boundaries, and to increment junc_hits.
    let left_ins = gexon.left_ins.max(0) as u32;
    let left_gap = gexon.left_gap.max(0) as u32;
    let right_ins = gexon.right_ins.max(0) as u32;
    let right_gap = gexon.right_gap.max(0) as u32;

    let cigar = match_info.align.cigar.as_mut().unwrap();

    // Handle start boundary
    if left_ins > 0 {
        if matches!(segment.status, ExonStatus::FirstExon | ExonStatus::OnlyExon) {
            if !has_left_clip {
                cigar.add_operation(left_ins, CigarOp::SoftClip);
                match_info.total_operations += left_ins as f64;
                match_info.prev_op = CigarOp::SoftClip;
            }
        } else if matches!(
            segment.status,
            ExonStatus::MiddleExon | ExonStatus::LastExon
        ) || has_left_clip
        {
            cigar.add_operation(left_ins, CigarOp::Ins);
            match_info.total_operations += left_ins as f64;
            if matches!(segment.status, ExonStatus::MiddleExon | ExonStatus::LastExon) {
                match_info.junc_misses += 1;
            }
            if match_info.prev_op == CigarOp::Del {
                match_info.total_coverage += left_ins as f64;
            } else if match_info.prev_op == CigarOp::Ins {
                match_info.total_operations += match_info.total_operations * 0.2;
            }
            match_info.prev_op = CigarOp::Ins;
        }
    } else if left_gap > 0 {
        if !first_match
            && (matches!(
                segment.status,
                ExonStatus::MiddleExon | ExonStatus::LastExon
            ) || has_left_clip)
        {
            cigar.add_operation(left_gap, CigarOp::Del);
            match_info.total_operations += left_gap as f64;
            match_info.ref_consumed += left_gap as i32;
            if matches!(segment.status, ExonStatus::MiddleExon | ExonStatus::LastExon) {
                match_info.junc_misses += 1;
            }
            if match_info.prev_op == CigarOp::Ins {
                match_info.total_coverage += left_gap as f64;
            } else if match_info.prev_op == CigarOp::Del {
                match_info.total_operations += match_info.total_operations * 0.2;
            }
            match_info.prev_op = CigarOp::Del;
        }
    } else {
        match_info.junc_hits += 1;
    }

    let overlap_start = qstart.max(gstart);
    let overlap_end = qend.min(gend);
    if overlap_end >= overlap_start {
        let match_length = overlap_end - overlap_start;
        cigar.add_operation(match_length, CigarOp::Match);
        match_info.total_operations += match_length as f64;
        match_info.total_coverage += match_length as f64;
        match_info.ref_consumed += match_length as i32;
        match_info.prev_op = CigarOp::Match;
    }

    // Handle end boundary
    if right_ins > 0 {
        if matches!(segment.status, ExonStatus::LastExon | ExonStatus::OnlyExon) {
            if !has_right_clip {
                cigar.add_operation(right_ins, CigarOp::SoftClip);
                match_info.total_operations += right_ins as f64;
                match_info.prev_op = CigarOp::SoftClip;
            }
        } else if matches!(
            segment.status,
            ExonStatus::FirstExon | ExonStatus::MiddleExon
        ) || has_right_clip
        {
            cigar.add_operation(right_ins, CigarOp::Ins);
            match_info.total_operations += right_ins as f64;
            if matches!(segment.status, ExonStatus::FirstExon | ExonStatus::MiddleExon) {
                match_info.junc_misses += 1;
            }
            if match_info.prev_op == CigarOp::Del {
                match_info.total_coverage += right_ins as f64;
            }
            match_info.prev_op = CigarOp::Ins;
        }
    } else if right_gap > 0 {
        if !last_match
            && (matches!(
                segment.status,
                ExonStatus::FirstExon | ExonStatus::MiddleExon
            ) || has_right_clip)
        {
            cigar.add_operation(right_gap, CigarOp::Del);
            match_info.total_operations += right_gap as f64;
            match_info.ref_consumed += right_gap as i32;
            if matches!(segment.status, ExonStatus::FirstExon | ExonStatus::MiddleExon) {
                match_info.junc_misses += 1;
            }
            if match_info.prev_op == CigarOp::Ins {
                match_info.total_coverage += right_gap as f64;
            }
            match_info.prev_op = CigarOp::Del;
        }
    } else {
        match_info.junc_hits += 1;
    }
}

pub fn build_cigar_ins(
    segment: &EvalSegment,
    k: usize,
    n: usize,
    match_info: &mut ExonChainMatch,
) {
    let qstart = segment.qexon.start;
    let qend = segment.qexon.end;
    let ins_length = qend.saturating_sub(qstart);
    let cigar = match_info.align.cigar.as_mut().unwrap();

    if k == 0 || k == (n - 1) {
        cigar.add_operation(ins_length, CigarOp::SoftClip);
        match_info.prev_op = CigarOp::SoftClip;
    } else {
        cigar.add_operation(ins_length, CigarOp::Ins);
        match_info.prev_op = CigarOp::Ins;
    }
    match_info.total_operations += ins_length as f64;
    match_info.total_coverage += ins_length as f64;
}

pub fn build_cigar_gap(segment: &EvalSegment, match_info: &mut ExonChainMatch) {
    let gexon = segment.gexon.as_ref().unwrap();
    let gap_length = gexon.end.saturating_sub(gexon.start);
    let cigar = match_info.align.cigar.as_mut().unwrap();
    cigar.add_operation(gap_length, CigarOp::Del);
    match_info.prev_op = CigarOp::Del;
    match_info.total_operations += gap_length as f64;
    match_info.total_coverage += gap_length as f64;
    match_info.ref_consumed += gap_length as i32;
}

pub fn build_cigar_clip(segment: &EvalSegment, match_info: &mut ExonChainMatch) {
    let cigar = match_info.align.cigar.as_mut().unwrap();
    for (len, op) in segment.cigar.ops.iter().copied() {
        cigar.add_operation(len, op);
        if matches!(op, CigarOp::MatchOverride | CigarOp::DelOverride) {
            match_info.ref_consumed += len as i32;
        }
    }
    match_info.align.clip_score += segment.score;
}

fn deterministic_index(name: &[u8], n: usize) -> usize {
    if n == 0 {
        return 0;
    }
    let state = RandomState::with_seeds(0, 0, 0, 0);
    (state.hash_one(name) % (n as u64)) as usize
}

pub fn filter_by_similarity(
    matches: &mut HashMap<Tid, ExonChainMatch>,
    config: &ReadEvaluationConfig,
    name: &[u8],
) {
    // C++ gates this whole similarity retain on `config.filter_by_similarity`
    // (= `similarity_threshold < 1.0`). Short reads (threshold 1.0) DISABLE it:
    // no match is dropped by similarity and `similarity_score` keeps its default,
    // so junc_hits never re-weights short-read placements (matching C++, which
    // also only rewrites the AS tag for long reads). When enabled (long reads),
    // drop sub-threshold matches and set `similarity_score = x²·(junc_hits+1)`.
    if config.filter_by_similarity {
        matches.retain(|_tid, m| {
            let similarity = if m.total_operations > 0.0 {
                m.total_coverage / m.total_operations
            } else {
                0.0
            };

            if similarity > config.similarity_threshold as f64 {
                let x = (similarity - config.similarity_threshold as f64)
                    / (1.0 - config.similarity_threshold as f64);
                // Discount transcripts the read's junctions disagree with: each
                // internal junction mismatch (a tolerated gap/ins) multiplies the
                // score by `junc_miss_discount` (<= 1), sharpening isoform choice.
                let miss_factor = if config.junc_miss_discount < 1.0 && m.junc_misses > 0 {
                    config.junc_miss_discount.powi(m.junc_misses)
                } else {
                    1.0
                };
                m.align.similarity_score = x * x * ((m.junc_hits + 1) as f64) * miss_factor;
                true
            } else {
                false
            }
        });
    }

    matches.retain(|_tid, m| {
        let pos = if m.align.strand == '+' {
            m.align.fwpos
        } else {
            m.align.rcpos
        } as i32;
        pos + m.ref_consumed <= m.transcript_len
    });

    if matches.is_empty() {
        return;
    }

    let mut best_tid: Option<Tid> = None;
    let mut best_score = f64::NEG_INFINITY;
    let mut count_at_best = 0;

    let mut tids: Vec<Tid> = matches.keys().copied().collect();
    tids.sort_unstable();

    for (hit_index, tid) in tids.iter().copied().enumerate() {
        let m = matches.get_mut(&tid).expect("tid must exist");
        m.align.hit_index = (hit_index + 1) as i32;
        if m.align.similarity_score > best_score {
            best_score = m.align.similarity_score;
            best_tid = Some(tid);
            count_at_best = 1;
        } else if m.align.similarity_score == best_score {
            count_at_best += 1;
        }
    }

    if let Some(tid) = best_tid
        && count_at_best == 1
        && let Some(m) = matches.get_mut(&tid)
    {
        m.align.primary_alignment = true;
        return;
    }

    // Tie or no unique best
    let idx = deterministic_index(name, tids.len());
    if let Some(tid) = tids.get(idx)
        && let Some(m) = matches.get_mut(tid)
    {
        m.align.primary_alignment = true;
    }
}

pub fn evaluate_exon_chains(
    read: &ReadAln,
    _id: ReadId,
    g2t: &G2TTree,
    config: ReadEvaluationConfig,
    long_reads: bool,
    shared_seq: Option<&[u8]>,
    ctx: &mut EvalContext,
) {
    ctx.matches.clear();
    let exon_count = read.segs.len();
    let refid = read.refid;

    let (has_left_clip, has_right_clip, n_left_clip, n_right_clip) = if long_reads {
        get_clips(read, &config).unwrap_or((false, false, 0, 0))
    } else {
        (false, false, 0, 0)
    };

    let strands_to_check = get_strands_to_check(read, long_reads);
    for strand in strands_to_check {
        ctx.data.clear();
        let mut failure = false;
        let mut has_gap = false;
        let mut start_ignore = EvalSegment::default();
        let mut end_ignore = EvalSegment::default();

        for j in 0..exon_count {
            get_intervals(
                ctx,
                read,
                j,
                exon_count,
                &config,
                g2t,
                refid,
                strand,
                has_left_clip,
                has_right_clip,
                &mut start_ignore,
                &mut end_ignore,
                &mut has_gap,
                &mut failure,
            );
            if failure {
                break;
            }
        }
        if failure {
            continue;
        }

        ctx.tids_buf.clear();
        ctx.tids_buf.extend(ctx.data.keys().copied());
        for tid in &ctx.tids_buf {
            let td = ctx.data.get_mut(tid).unwrap();
            if td.elim {
                continue;
            }
            correct_for_gaps(td, *tid, &config, g2t, strand, refid, long_reads);
        }

        if long_reads && config.use_fasta {
            for tid in &ctx.tids_buf {
                let td = ctx.data.get_mut(tid).unwrap();
                if td.elim {
                    continue;
                }
                if td.has_left_clip {
                    if n_left_clip >= 5 {
                        left_clip_rescue(
                            td,
                            strand,
                            g2t,
                            refid,
                            *tid,
                            n_left_clip,
                            &config,
                            read,
                            shared_seq,
                            &mut ctx.aligner,
                            &mut ctx.sw_bufs,
                        );
                    } else {
                        td.has_left_clip = false;
                    }
                }

                if td.has_right_clip {
                    if n_right_clip >= 5 {
                        right_clip_rescue(
                            td,
                            strand,
                            g2t,
                            refid,
                            *tid,
                            n_right_clip,
                            &config,
                            read,
                            shared_seq,
                            &mut ctx.aligner,
                            &mut ctx.sw_bufs,
                        );
                    } else {
                        td.has_right_clip = false;
                    }
                }
            }
        }

        for (tid, td) in ctx.data.iter_mut() {
            if td.elim {
                continue;
            }

            let n_segments = td.segments.len();
            let mut n_ins_exon: u32 = 0;
            let mut match_created = false;
            let mut first_match_idx: i32 = -1;
            let mut last_match_idx: i32 = -1;
            let mut prev_s = 0;
            let mut prev_e = 0;

            for segment in td.segments.iter() {
                if !match_created && segment.status == ExonStatus::InsExon {
                    n_ins_exon += segment.qexon.end - segment.qexon.start;
                    // C++ skips InsExon entirely (has_gexon=false) — do NOT
                    // increment first_match_idx/last_match_idx here.
                    continue;
                }

                if let Some(gexon) = segment.gexon.as_ref() {
                    if gexon.start == prev_s && gexon.end == prev_e {
                        td.elim = true;
                        break;
                    }
                    prev_s = gexon.start;
                    prev_e = gexon.end;
                }

                // Match C++: use has_gexon (not status check) for both conditions.
                // Clip segments have has_gexon=false (matching C++), so they're
                // skipped here just like in C++.
                if !match_created && segment.has_gexon {
                    let gexon = segment.gexon.as_ref().unwrap();
                    let read_strand_for_match = if long_reads
                        && read.segs.len() == 1
                        && !read.strand_from_tag
                    {
                        strand
                    } else {
                        read.strand
                    };
                    td.match_info =
                        create_match(gexon, strand, read_strand_for_match, n_ins_exon);
                    match_created = true;
                    first_match_idx += 1;
                    last_match_idx += 1;
                } else if match_created
                    && segment.has_gexon
                    && segment.status != ExonStatus::InsExon
                {
                    last_match_idx += 1;
                    if strand == '-'
                        && let Some(gexon) = segment.gexon.as_ref()
                    {
                        td.match_info.align.rcpos = gexon.pos;
                    }
                }
            }

            if td.elim {
                continue;
            }

            if td.match_info.align.cigar.is_none() {
                td.elim = true;
                continue;
            }

            let mut k_build: i32 = -1;
            for segment in td.segments.iter() {
                k_build += 1;
                let first_match = k_build == first_match_idx;
                let last_match = k_build == last_match_idx;

                let is_match = matches!(
                    segment.status,
                    ExonStatus::FirstExon
                        | ExonStatus::MiddleExon
                        | ExonStatus::LastExon
                        | ExonStatus::OnlyExon
                );
                let is_ins = segment.status == ExonStatus::InsExon;
                let is_gap = segment.status == ExonStatus::GapExon;
                let is_clip = matches!(
                    segment.status,
                    ExonStatus::LeftClipExon | ExonStatus::RightClipExon
                );

                if is_match {
                    let has_left_clip = td.has_left_clip;
                    let has_right_clip = td.has_right_clip;
                    build_cigar_match(
                        segment,
                        has_left_clip,
                        has_right_clip,
                        &mut td.match_info,
                        first_match,
                        last_match,
                    );
                } else if is_ins {
                    build_cigar_ins(segment, k_build as usize, n_segments, &mut td.match_info);
                    if k_build == 0 || k_build == (n_segments as i32 - 1) {
                        td.match_info.junc_hits -= 1;
                    } else {
                        td.match_info.junc_hits -= 2;
                    }
                } else if is_gap {
                    build_cigar_gap(segment, &mut td.match_info);
                    td.match_info.junc_hits -= 2;
                } else if is_clip {
                    build_cigar_clip(segment, &mut td.match_info);
                }
            }

            if td.match_info.junc_hits < 0 {
                td.match_info.junc_hits = 0;
            }

            if !td.elim {
                ctx.matches.insert(*tid, std::mem::take(&mut td.match_info));
            }
        }
    }

    if !ctx.matches.is_empty() {
        filter_by_similarity(&mut ctx.matches, &config, &read.name);
    }
}
#[derive(Debug, Default)]
pub struct ReadEvaluator {
    pub lr: bool,
    pub lr_hq: bool,
    pub strict: bool,
    pub use_fasta: bool,
    pub max_clip: Option<u8>,
    pub max_ins: Option<u8>,
    pub max_junc_gap: Option<u8>,
    pub similarity_threshold: Option<f32>,
    pub small_exon_size: Option<u8>,
    /// Per-junction-mismatch similarity discount in `(0, 1]`; `None` => 1.0 (off).
    pub junc_miss_discount: Option<f64>,
}

impl ReadEvaluator {
    fn build_config(&self) -> ReadEvaluationConfig {
        let (default_max_clip, default_max_ins, default_max_junc_gap, default_similarity_threshold, default_small_exon_size) = 
            match (self.lr, self.lr_hq, self.strict) {
                (_, _, true) => (0,   0,  0, 0.99_f32,  0),
                (_, true, _) => (5,  10, 10, 0.90_f32,  0),
                (true, _, _) => (40, 40, 40, 0.60_f32, 35),
                // Short-read default matches C++ (SIM_THR.value_or(1.0)): threshold
                // 1.0 is a SENTINEL that DISABLES the similarity filter
                // (`filter_by_similarity = threshold < 1.0`), NOT a "require a
                // perfect match" gate. With the filter gated off (below), reads are
                // no longer dropped by the `similarity > 1.0` test — that gate is
                // never evaluated — and short-read placements are not re-weighted by
                // junc_hits, matching C++ (which only rewrites AS for long reads).
                _            => (5,   0,  0, 1.0_f32,  0),
            };

        let small_exon_size = self.small_exon_size.unwrap_or(default_small_exon_size);
        let similarity_threshold = self.similarity_threshold.unwrap_or(default_similarity_threshold);

        ReadEvaluationConfig {
            max_clip:               self.max_clip.unwrap_or(default_max_clip) as u32,
            max_ins:                self.max_ins.unwrap_or(default_max_ins) as u32,
            max_junc_gap:           self.max_junc_gap.unwrap_or(default_max_junc_gap) as u32,
            similarity_threshold,
            // C++: filter_by_similarity = (similarity_threshold < 1.0).
            filter_by_similarity:   similarity_threshold < 1.0,
            ignore_small_exons:     small_exon_size > 0,
            small_exon_size:        small_exon_size as u32,
            print:                  false,
            name:                   Vec::new(),
            soft_clips:             true,
            strict:                 self.strict,
            use_fasta:              self.use_fasta,
            junc_miss_discount:     self.junc_miss_discount.unwrap_or(1.0),
        }
    }

    pub fn evaluate(
        &self,
        read: &ReadAln,
        id: ReadId,
        g2t: &G2TTree,
        seq: Option<&[u8]>,
        ctx: &mut EvalContext,
    ) {
        let config = self.build_config();
        evaluate_exon_chains(read, id, g2t, config, self.lr || self.lr_hq, seq, ctx);
    }
}

#[cfg(test)]
mod tests {
    use super::{Cigar, CigarOp, ReadEvaluator};

    #[test]
    fn short_read_similarity_filter_is_disabled() {
        // C++ parity: a default (short-read) ReadEvaluator resolves to threshold
        // 1.0, which DISABLES the similarity filter (`filter_by_similarity =
        // threshold < 1.0 = false`). The filter being off is what prevents the
        // `similarity > threshold` gate from ever running (so short reads are NOT
        // dropped — the earlier "projected 0" bug) AND leaves similarity_score
        // unweighted (so short-read placements aren't re-scored by junc_hits,
        // matching C++ which only rewrites AS for long reads).
        let short = ReadEvaluator::default().build_config();
        assert_eq!(short.similarity_threshold, 1.0, "short-read threshold is 1.0");
        assert!(!short.filter_by_similarity, "short-read similarity filter is disabled");

        // Long reads keep the filter enabled (threshold 0.60 < 1.0).
        let long = ReadEvaluator { lr: true, ..Default::default() }.build_config();
        assert!(long.filter_by_similarity, "long-read similarity filter is enabled");
        assert!((long.similarity_threshold - 0.60).abs() < 1e-6);
    }

    #[test]
    fn query_consumed_counts_only_query_ops() {
        // 10M 3I 5D 4S 2N: query-consuming = M(10) + I(3) = 13.
        let mut c = Cigar::default();
        c.ops.push((10, CigarOp::Match));
        c.ops.push((3, CigarOp::Ins));
        c.ops.push((5, CigarOp::Del));
        c.ops.push((4, CigarOp::SoftClip));
        c.ops.push((2, CigarOp::RefSkip));
        assert_eq!(c.query_consumed(), 13);
    }

    #[test]
    fn query_consumed_counts_equal_diff_and_overrides() {
        // 6= 1X 2/ (InsOverride) 3, (MatchOverride): all consume query = 12.
        let mut c = Cigar::default();
        c.ops.push((6, CigarOp::Equal));
        c.ops.push((1, CigarOp::Diff));
        c.ops.push((2, CigarOp::InsOverride));
        c.ops.push((3, CigarOp::MatchOverride));
        assert_eq!(c.query_consumed(), 12);
    }
}
