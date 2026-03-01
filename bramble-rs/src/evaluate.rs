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
    pub max_gap: u32,
    pub ignore_small_exons: bool,
    pub small_exon_size: u32,
    pub max_junc_gap: u32,
    pub similarity_threshold: f32,
    #[allow(dead_code)]
    pub print: bool,
    pub name: String,
    pub soft_clips: bool,
    pub strict: bool,
    pub use_fasta: bool,
}

impl ReadEvaluationConfig {
    /// Parameters matching C++ ShortReadEvaluator.
    pub fn short_read() -> Self {
        Self {
            max_clip: 5,
            max_ins: 0,
            max_gap: 0,
            ignore_small_exons: false,
            small_exon_size: 0,
            max_junc_gap: 0,
            similarity_threshold: 0.90,
            print: false,
            name: String::new(),
            soft_clips: true,
            strict: false,
            use_fasta: false,
        }
    }

    /// Parameters matching C++ LongReadEvaluator.
    #[allow(dead_code)]
    pub fn long_read() -> Self {
        Self {
            max_clip: 40,
            max_ins: 40,
            max_gap: 40,
            ignore_small_exons: true,
            small_exon_size: 35,
            max_junc_gap: 40,
            similarity_threshold: 0.60,
            print: false,
            name: String::new(),
            soft_clips: true,
            strict: false,
            use_fasta: false,
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
    pub has_left_ins: bool,
    pub has_right_ins: bool,
    pub n_left_ins: i32,
    pub n_right_ins: i32,
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
    pub name: String,
}

/// Reusable scratch space for one evaluation pass.
/// Allocate once per worker thread; cleared between reads automatically.
pub(crate) struct EvalContext {
    /// Final output: transcript → chain match.  Filled by evaluate(); drained by caller.
    pub(crate) matches: HashMap<Tid, ExonChainMatch>,
    /// Per-strand working map: tid → accumulated TidData.
    pub(crate) data: HashMap<Tid, TidData>,
    /// Working set of matched tids for the current exon.
    pub(crate) candidate_tids: HashSet<Tid>,
    /// Scratch buffer for guide exon queries.
    pub(crate) guide_exons: HashMap<Tid, GuideExon>,
}

impl EvalContext {
    pub(crate) fn new() -> Self {
        Self {
            matches: HashMap::new(),
            data: HashMap::new(),
            candidate_tids: HashSet::new(),
            guide_exons: HashMap::new(),
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

    let ops: Vec<_> = read.cigar_ops.to_vec();

    if ops.is_empty() {
        return Ok((false, false, 0, 0));
    }

    let (left0_len, left0_kind) = ops.first().unwrap();
    let left1 = ops.get(1);
    if *left0_kind == CigarKind::HardClip {
        if let Some((op_len, op_kind)) = left1
            && *op_kind == CigarKind::SoftClip
        {
            has_left_clip = config.soft_clips;
            n_left_clip = *op_len;
        }
    } else if *left0_kind == CigarKind::SoftClip {
        has_left_clip = config.soft_clips;
        n_left_clip = *left0_len;
    }

    let (right0_len, right0_kind) = ops.last().unwrap();
    let right1 = if ops.len() >= 2 { ops.get(ops.len() - 2) } else { None };
    if *right0_kind == CigarKind::HardClip {
        if let Some((op_len, op_kind)) = right1
            && *op_kind == CigarKind::SoftClip
        {
            has_right_clip = config.soft_clips;
            n_right_clip = *op_len;
        }
    } else if *right0_kind == CigarKind::SoftClip {
        has_right_clip = config.soft_clips;
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
                    has_left_ins: gexon.left_ins > 0,
                    n_left_ins: gexon.left_ins,
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
                && let Some(tid_data) = ctx.data.get_mut(&tid)
            {
                tid_data.has_right_ins = gexon.right_ins > 0;
                tid_data.n_right_ins = gexon.right_ins;
                if gexon.right_gap > 0 {
                    *has_gap = true;
                }
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
    start_ignore: &EvalSegment,
    strand: char,
    g2t: &G2TTree,
    refid: RefId,
    tid: Tid,
    n_left_clip: u32,
    n_left_ins: i32,
    config: &ReadEvaluationConfig,
    read: &ReadAln,
    shared_seq: Option<&[u8]>,
) {
    if !config.use_fasta {
        return;
    }
    if start_ignore.is_valid {
        tid_data.has_left_clip = false;
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

    let total_bases = (n_left_clip as i32 + n_left_ins).max(0) as usize;
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
    // NOTE: structured as a labeled block rather than a loop so that early exits via
    // `break 'left_clip` are explicit.  The loop was intentionally single-iteration;
    // the label preserves the ability to extend to multi-exon rescue in the future.
    'left_clip: {
        if remaining_qseq.is_empty() {
            break 'left_clip;
        }

        let (left_gexon, next_gseq) = if strand == '+' {
            if !current_gexon.has_prev {
                tid_data.has_left_clip = false;
                return;
            }
            let left = g2t.get_guide_exon_for_tid(
                refid,
                strand,
                tid,
                current_gexon.prev_start,
                current_gexon.prev_end,
            );
            let Some(left) = left else {
                tid_data.has_left_clip = false;
                return;
            };
            let next = left.seq.clone();
            (left, next)
        } else {
            if !current_gexon.has_next {
                tid_data.has_left_clip = false;
                return;
            }
            let left = g2t.get_guide_exon_for_tid(
                refid,
                strand,
                tid,
                current_gexon.next_start,
                current_gexon.next_end,
            );
            let Some(left) = left else {
                tid_data.has_left_clip = false;
                return;
            };
            let next = left.seq.clone();
            (left, next)
        };

        if next_gseq.is_empty() {
            tid_data.has_left_clip = false;
            return;
        }
        let result = sw::smith_waterman(&remaining_qseq, next_gseq.as_bytes(), Anchor::End);
        if result.accuracy <= 0.70 || result.score < 10 {
            tid_data.has_left_clip = false;
            return;
        }

        let qstart = left_gexon.start + result.start_j as u32;
        let qend = left_gexon.start + result.end_j as u32;
        let mut left_gexon = left_gexon;
        left_gexon.pos = (result.pos as u32) + left_gexon.pos_start;

        if n_left_ins > 0 {
            segment.qexon.start = gexon.start;
        }

        let left_clip = EvalSegment {
            is_valid: true,
            has_qexon: true,
            has_gexon: true,
            gexon: Some(left_gexon.clone()),
            qexon: alignment::Segment { start: qstart, end: qend },
            status: ExonStatus::LeftClipExon,
            is_small_exon: remaining_qseq.len() as u32 <= config.small_exon_size,
            cigar: result.cigar.clone(),
            score: result.score,
        };
        tid_data.segments.insert(0, left_clip);
        tid_data.has_left_clip = true;

        if result.start_i > 1 {
            remaining_qseq.truncate(result.start_i as usize);
        } else {
            remaining_qseq.clear();
        }
    }

    if !remaining_qseq.is_empty()
        && let Some(first) = tid_data.segments.first_mut()
    {
        first
            .cigar
            .prepend_operation(remaining_qseq.len() as u32, CigarOp::ClipOverride);
    }
}

#[allow(clippy::too_many_arguments)]
pub fn right_clip_rescue(
    tid_data: &mut TidData,
    end_ignore: &EvalSegment,
    strand: char,
    g2t: &G2TTree,
    refid: RefId,
    tid: Tid,
    n_right_clip: u32,
    n_right_ins: i32,
    config: &ReadEvaluationConfig,
    read: &ReadAln,
    shared_seq: Option<&[u8]>,
) {
    if !config.use_fasta {
        return;
    }
    if end_ignore.is_valid {
        tid_data.has_right_clip = false;
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

    let total_bases = (n_right_clip as i32 + n_right_ins).max(0) as usize;
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
    // NOTE: structured as a labeled block — see left_clip_rescue for rationale.
    'right_clip: {
        if remaining_qseq.is_empty() {
            break 'right_clip;
        }

        let right_gexon = if strand == '+' {
            if !current_gexon.has_next {
                tid_data.has_right_clip = false;
                return;
            }
            g2t.get_guide_exon_for_tid(
                refid,
                strand,
                tid,
                current_gexon.next_start,
                current_gexon.next_end,
            )
        } else {
            if !current_gexon.has_prev {
                tid_data.has_right_clip = false;
                return;
            }
            g2t.get_guide_exon_for_tid(
                refid,
                strand,
                tid,
                current_gexon.prev_start,
                current_gexon.prev_end,
            )
        };

        let Some(mut right_gexon) = right_gexon else {
            tid_data.has_right_clip = false;
            return;
        };
        if right_gexon.seq.is_empty() {
            tid_data.has_right_clip = false;
            return;
        }

        let result = sw::smith_waterman(&remaining_qseq, right_gexon.seq.as_bytes(), Anchor::Start);
        if result.accuracy <= 0.70 || result.score < 10 {
            tid_data.has_right_clip = false;
            return;
        }

        let qstart = right_gexon.start + result.start_j as u32;
        let qend = right_gexon.start + result.end_j as u32;
        right_gexon.pos = (result.pos as u32) + right_gexon.pos_start;

        if n_right_ins > 0 {
            segment.qexon.end = gexon.end;
        }

        let right_clip = EvalSegment {
            is_valid: true,
            has_qexon: true,
            has_gexon: true,
            gexon: Some(right_gexon.clone()),
            qexon: alignment::Segment { start: qstart, end: qend },
            status: ExonStatus::RightClipExon,
            is_small_exon: remaining_qseq.len() as u32 <= config.small_exon_size,
            cigar: result.cigar.clone(),
            score: result.score,
        };
        tid_data.segments.push(right_clip);
        tid_data.has_right_clip = true;

        let query_consumed = result.end_i + 1;
        if query_consumed < remaining_qseq.len() as i32 {
            remaining_qseq.drain(..query_consumed as usize);
        } else {
            remaining_qseq.clear();
        }
    }

    if !remaining_qseq.is_empty()
        && let Some(last) = tid_data.segments.last_mut()
    {
        last
            .cigar
            .add_operation(remaining_qseq.len() as u32, CigarOp::ClipOverride);
    }
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
    let mut out: Vec<EvalSegment> = Vec::with_capacity(tid_data.segments.len() + 4);
    let mut prev_gexon: Option<GuideExon> = None;

    for segment in tid_data.segments.iter() {
        if tid_data.elim {
            return;
        }
        if !segment.has_gexon {
            out.push(segment.clone());
            continue;
        }

        let gexon = segment.gexon.clone().unwrap();

        if let Some(prev) = prev_gexon.clone() {
            let gap = gexon.exon_id.saturating_sub(prev.exon_id);

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
                        || gap_end.saturating_sub(gap_start) > config.max_gap
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

        prev_gexon = Some(gexon.clone());
        out.push(segment.clone());
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

    let cigar = match_info.align.cigar.as_mut().unwrap();

    // Handle start boundary
    if qstart < gstart {
        let overhang = gstart - qstart;
        if matches!(segment.status, ExonStatus::FirstExon | ExonStatus::OnlyExon) {
            if !has_left_clip {
                cigar.add_operation(overhang, CigarOp::SoftClip);
                match_info.total_operations += overhang as f64;
                match_info.prev_op = CigarOp::SoftClip;
            }
        } else if matches!(
            segment.status,
            ExonStatus::MiddleExon | ExonStatus::LastExon
        ) || has_left_clip
        {
            cigar.add_operation(overhang, CigarOp::Ins);
            match_info.total_operations += overhang as f64;
            if match_info.prev_op == CigarOp::Del {
                match_info.total_coverage += overhang as f64;
            } else if match_info.prev_op == CigarOp::Ins {
                match_info.total_operations += match_info.total_operations * 0.2;
            }
            match_info.prev_op = CigarOp::Ins;
        }
    } else if gstart < qstart {
        if !first_match
            && (matches!(
                segment.status,
                ExonStatus::MiddleExon | ExonStatus::LastExon
            ) || has_left_clip)
        {
            let overhang = qstart - gstart;
            cigar.add_operation(overhang, CigarOp::Del);
            match_info.total_operations += overhang as f64;
            match_info.ref_consumed += overhang as i32;
            if match_info.prev_op == CigarOp::Ins {
                match_info.total_coverage += overhang as f64;
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
    if gend < qend {
        let overhang = qend - gend;
        if matches!(segment.status, ExonStatus::LastExon | ExonStatus::OnlyExon) {
            if !has_right_clip {
                cigar.add_operation(overhang, CigarOp::SoftClip);
                match_info.total_operations += overhang as f64;
                match_info.prev_op = CigarOp::SoftClip;
            }
        } else if matches!(
            segment.status,
            ExonStatus::FirstExon | ExonStatus::MiddleExon
        ) || has_right_clip
        {
            cigar.add_operation(overhang, CigarOp::Ins);
            match_info.total_operations += overhang as f64;
            if match_info.prev_op == CigarOp::Del {
                match_info.total_coverage += overhang as f64;
            }
            match_info.prev_op = CigarOp::Ins;
        }
    } else if qend < gend {
        if !last_match
            && (matches!(
                segment.status,
                ExonStatus::FirstExon | ExonStatus::MiddleExon
            ) || has_right_clip)
        {
            let overhang = gend - qend;
            cigar.add_operation(overhang, CigarOp::Del);
            match_info.total_operations += overhang as f64;
            match_info.ref_consumed += overhang as i32;
            if match_info.prev_op == CigarOp::Ins {
                match_info.total_coverage += overhang as f64;
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

fn deterministic_index(name: &str, n: usize) -> usize {
    if n == 0 {
        return 0;
    }
    let state = RandomState::with_seeds(0, 0, 0, 0);
    (state.hash_one(name) % (n as u64)) as usize
}

pub fn filter_by_similarity(
    matches: &mut HashMap<Tid, ExonChainMatch>,
    config: &ReadEvaluationConfig,
) {
    matches.retain(|_tid, m| {
        let similarity = if m.total_operations > 0.0 {
            m.total_coverage / m.total_operations
        } else {
            0.0
        };

        if similarity > config.similarity_threshold as f64 {
            let x = (similarity - config.similarity_threshold as f64)
                / (1.0 - config.similarity_threshold as f64);
            m.align.similarity_score = x * x * ((m.junc_hits + 1) as f64);
            true
        } else {
            false
        }
    });

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
    let idx = deterministic_index(&config.name, tids.len());
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
    mut config: ReadEvaluationConfig,
    long_reads: bool,
    shared_seq: Option<&[u8]>,
    ctx: &mut EvalContext,
) {
    ctx.matches.clear();
    let exon_count = read.segs.len();
    let refid = read.refid;
    config.name = read.name.clone();

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

        let tids: Vec<Tid> = ctx.data.keys().copied().collect();
        for tid in &tids {
            let td = ctx.data.get_mut(tid).unwrap();
            if td.elim {
                continue;
            }
            correct_for_gaps(td, *tid, &config, g2t, strand, refid, long_reads);
        }

        if long_reads && config.use_fasta {
            for tid in &tids {
                let td = ctx.data.get_mut(tid).unwrap();
                if td.elim {
                    continue;
                }
                if td.has_left_ins {
                    if n_left_clip >= 5 {
                        left_clip_rescue(
                            td,
                            &start_ignore,
                            strand,
                            g2t,
                            refid,
                            *tid,
                            n_left_clip,
                            td.n_left_ins,
                            &config,
                            read,
                            shared_seq,
                        );
                    } else {
                        td.has_left_clip = false;
                    }
                } else if td.has_left_clip {
                    if n_left_clip >= 5 {
                        left_clip_rescue(
                            td,
                            &start_ignore,
                            strand,
                            g2t,
                            refid,
                            *tid,
                            n_left_clip,
                            0,
                            &config,
                            read,
                            shared_seq,
                        );
                    } else {
                        td.has_left_clip = false;
                    }
                }

                if td.has_right_ins {
                    if n_right_clip >= 5 {
                        right_clip_rescue(
                            td,
                            &end_ignore,
                            strand,
                            g2t,
                            refid,
                            *tid,
                            n_right_clip,
                            td.n_right_ins,
                            &config,
                            read,
                            shared_seq,
                        );
                    } else {
                        td.has_right_clip = false;
                    }
                } else if td.has_right_clip {
                    if n_right_clip >= 5 {
                        right_clip_rescue(
                            td,
                            &end_ignore,
                            strand,
                            g2t,
                            refid,
                            *tid,
                            n_right_clip,
                            0,
                            &config,
                            read,
                            shared_seq,
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
                    first_match_idx += 1;
                    last_match_idx += 1;
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

                if !match_created && segment.status != ExonStatus::InsExon {
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
                } else if match_created && segment.status != ExonStatus::InsExon {
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
                ctx.matches.insert(*tid, td.match_info.clone());
            }
        }
    }

    if !ctx.matches.is_empty() {
        filter_by_similarity(&mut ctx.matches, &config);
    }
}
#[derive(Debug, Default)]
pub struct ReadEvaluator {
    pub long_reads: bool,
    pub use_fasta: bool,
}

impl ReadEvaluator {
    pub fn evaluate(
        &self,
        read: &ReadAln,
        id: ReadId,
        g2t: &G2TTree,
        seq: Option<&[u8]>,
        ctx: &mut EvalContext,
    ) {
        let (max_clip, max_ins, max_gap, max_junc_gap, similarity_threshold, ignore_small_exons, small_exon_size) =
            if self.long_reads {
                (40, 40, 40, 40, 0.60_f32, true, 35)
            } else {
                (5, 0, 0, 0, 0.90_f32, false, 0)
            };

        let config = ReadEvaluationConfig {
            max_clip,
            max_ins,
            max_gap,
            ignore_small_exons,
            small_exon_size,
            max_junc_gap,
            similarity_threshold,
            print: false,
            name: String::new(),
            soft_clips: true,
            strict: false,
            use_fasta: self.use_fasta,
        };

        evaluate_exon_chains(read, id, g2t, config, self.long_reads, seq, ctx);
    }
}
