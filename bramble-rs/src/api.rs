//! Public library API for projecting spliced genomic alignments into transcriptome space.
//!
//! The entry point is [`project_group`]: given a slice of [`GenomicAlignment`]s for the
//! same query name and a pre-built [`G2TTree`] index, it returns one
//! [`ProjectedAlignment`] per successful transcript match.
//!
//! # Quick start
//!
//! ```no_run
//! use bramble_rs::{GenomicAlignment, ProjectionConfig, project_group};
//! use bramble_rs::annotation::load_transcripts;
//! use bramble_rs::g2t::{G2TTree, build_g2t};
//! use std::collections::HashMap;
//!
//! // 1. Build the index once from a GTF/GFF annotation.
//! // let transcripts = load_transcripts("annotation.gtf")?;
//! // let refname_to_id: HashMap<String, usize> = /* from your BAM header */;
//! // let index = build_g2t(&transcripts, &refname_to_id, None)?;
//!
//! // 2. Choose short-read or long-read parameters.
//! // let config = ProjectionConfig { long_reads: false, use_fasta: false };
//!
//! // 3. For each read-name group, project all alignments at once.
//! // let alns: Vec<GenomicAlignment> = /* built from noodles RecordBuf or minimap2-rs Mapping */;
//! // let projected: Vec<ProjectedAlignment> = project_group(&alns, &index, &config);
//! ```
//!
//! # minimap2-rs integration sketch
//!
//! ```no_run
//! # use bramble_rs::GenomicAlignment;
//! // let mapping: minimap2::Mapping = /* ... */;
//! // let aln = GenomicAlignment {
//! //     query_name: read_name.clone(),
//! //     ref_id: mapping.target_name
//! //         .and_then(|n| refname_to_id.get(n).copied())
//! //         .unwrap_or(-1) as i32,
//! //     ref_start: mapping.target_start as i64 + 1, // convert 0-based → 1-based
//! //     is_reverse: mapping.strand == minimap2::Strand::Reverse,
//! //     cigar: mapping.alignment.as_ref()
//! //         .map(|a| a.cigar.iter().map(|&(len, op)| (len, op)).collect())
//! //         .unwrap_or_default(),
//! //     sequence: Some(read_sequence.clone()),
//! //     is_paired: false, is_first_in_pair: false,
//! //     xs_strand: None, ts_strand: mapping.strand_from_tag,
//! //     hit_index: 0, mate_ref_id: None, mate_ref_start: None,
//! //     mate_is_unmapped: false, read_len: read_sequence.len(),
//! // };
//! ```

use crate::evaluate::{ReadAln, ReadEvaluator, ExonChainMatch};
use crate::g2t::G2TTree;
use crate::pipeline::{
    ReadEval, OutputEntry,
    find_mate_pairs, assign_pair_order,
    build_paired_groups, build_unpaired_groups,
    assign_hit_indices, align_pos, compute_template_length,
    sam_op_to_kind, segs_from_ops,
};
use crate::types::{ReadId, Tid};
use noodles::sam::alignment::record::cigar::op::Kind as CigarKind;
use noodles::sam::alignment::record::Flags;
use std::collections::HashMap;

/// A genomic alignment ready for projection into transcriptome coordinates.
///
/// Construct one of these from your alignment source (noodles `RecordBuf`,
/// minimap2-rs `Mapping`, etc.) and pass a slice of them to [`project_group`].
///
/// # Coordinate conventions
///
/// - `ref_id` is a **0-based** index into the same reference-sequence table that
///   was used to build the [`G2TTree`] index.  It must match the integer the
///   index was built with — typically the position of the chromosome name in the
///   BAM/CRAM header or the minimap2 target list.
/// - `ref_start` is **1-based** (SAM convention).
#[derive(Debug, Clone)]
pub struct GenomicAlignment {
    /// Read/query name shared by all alignments in the group.
    pub query_name: String,
    /// 0-based reference sequence index (must match the index used to build [`G2TTree`]).
    pub ref_id: i32,
    /// 1-based alignment start position on the reference (SAM `POS`).
    pub ref_start: i64,
    /// `true` if the read maps to the reverse strand.
    pub is_reverse: bool,
    /// CIGAR as `(length, SAM op code)` pairs.
    ///
    /// SAM op codes: `0`=M, `1`=I, `2`=D, `3`=N (intron skip), `4`=S, `5`=H,
    /// `6`=P, `7`==, `8`=X.
    pub cigar: Vec<(u32, u8)>,
    /// Query sequence bytes (ASCII bases, e.g. `b"ACGT..."`).
    ///
    /// Required only when soft-clip rescue is enabled (`use_fasta = true` in
    /// [`ProjectionConfig`]).  May be `None` otherwise.
    pub sequence: Option<Vec<u8>>,
    /// `true` if the read comes from a paired-end library (SAM `SEGMENTED` flag).
    pub is_paired: bool,
    /// `true` if this is read 1 of a pair (SAM `FIRST_SEGMENT` flag).
    pub is_first_in_pair: bool,
    /// Library strand inferred from the SAM `XS` tag (`'+'` or `'-'`).
    ///
    /// Used for short-read strand inference.  Set to `None` if the tag is absent.
    pub xs_strand: Option<char>,
    /// Library strand inferred from the minimap2 `ts` tag (`'+'` or `'-'`).
    ///
    /// Used for long-read strand inference.  The value is the *template* strand
    /// (same convention as minimap2); bramble flips it when the read is reverse-
    /// complemented.  Set to `None` if the tag is absent.
    pub ts_strand: Option<char>,
    /// Value of the SAM `HI` tag from the input BAM (0 if the tag is absent).
    ///
    /// Used to match mate pairs when the same query has multiple secondary
    /// alignments.
    pub hit_index: i32,
    /// 0-based reference sequence index of the mate.
    ///
    /// `None` when the read is unpaired or the mate is unmapped.
    pub mate_ref_id: Option<i32>,
    /// 1-based alignment start of the mate (SAM `MPOS`).
    ///
    /// `None` when the read is unpaired or the mate is unmapped.
    pub mate_ref_start: Option<i64>,
    /// `true` if the mate is unmapped (SAM `MATE_UNMAPPED` flag).
    pub mate_is_unmapped: bool,
    /// Query read length in bases.
    ///
    /// Used to compute NH when [`sequence`](GenomicAlignment::sequence) is
    /// absent.  Set to 0 if unknown (the sequence length will be used instead).
    pub read_len: usize,
}

/// The result of projecting one alignment into transcriptome space.
///
/// Each value returned by [`project_group`] corresponds to one transcript that
/// the read successfully maps to.  Multiple `ProjectedAlignment`s may share the
/// same [`input_index`](ProjectedAlignment::input_index) when one genomic
/// alignment maps to several overlapping transcripts.
#[derive(Debug, Clone)]
pub struct ProjectedAlignment {
    /// Index of the matching transcript inside the [`G2TTree`].
    ///
    /// This is a dense 0-based integer assigned by [`build_g2t`](crate::g2t::build_g2t)
    /// in the order transcripts were inserted.
    pub transcript_id: u32,
    /// 1-based start position on the transcript (forward-strand coordinates).
    pub transcript_start: u32,
    /// `true` if the read is on the reverse strand of the transcript.
    pub is_reverse: bool,
    /// Alignment similarity score in `[0, 1]` (higher is better).
    ///
    /// Computed as the fraction of query bases that align to annotated exonic
    /// sequence without exceeding the clip/gap tolerances.
    pub similarity_score: f64,
    /// Total number of transcript hits for this read group (`NH` tag value).
    pub nh: u32,
    /// 1-based index of this particular hit (`HI` tag value).
    pub hi: u32,
    /// `true` if this is the primary (best-scoring) alignment.
    pub is_primary: bool,
    /// `true` if the mate read maps to the same transcript.
    pub same_transcript_as_mate: bool,
    /// Observed template length (`TLEN`); `0` when mates map to different
    /// transcripts or when the read is unpaired.
    pub insert_size: i32,
    /// Index of the input [`GenomicAlignment`] this result was derived from.
    ///
    /// Use this to correlate projected results back to the original BAM records
    /// or minimap2 mappings.
    pub input_index: usize,
}

/// Evaluation parameters for [`project_group`].
///
/// Use [`ProjectionConfig::short_read`] or [`ProjectionConfig::long_read`] for
/// pre-set configurations that match the C++ bramble defaults, or build a
/// custom config for advanced use.
#[derive(Debug, Clone)]
pub struct ProjectionConfig {
    /// Use long-read evaluation parameters (more permissive clip/gap tolerances
    /// and a lower similarity threshold).
    ///
    /// Set `true` for PacBio / Oxford Nanopore data, `false` for Illumina.
    pub long_reads: bool,
    /// Whether a genome FASTA sequence is available for soft-clip rescue
    /// alignment.
    ///
    /// When `true`, reads with soft-clipped bases are re-aligned against the
    /// reference sequence to recover additional exon coverage.  Requires that
    /// `sequence` is set on the [`GenomicAlignment`].
    pub use_fasta: bool,
}

impl ProjectionConfig {
    /// Short-read (Illumina) defaults — matches C++ `ShortReadEvaluator`.
    pub fn short_read() -> Self {
        Self { long_reads: false, use_fasta: false }
    }

    /// Long-read (PacBio / ONT) defaults — matches C++ `LongReadEvaluator`.
    pub fn long_read() -> Self {
        Self { long_reads: true, use_fasta: false }
    }
}

/// Project a group of alignments for the same query name into transcriptome coordinates.
///
/// # Parameters
///
/// * `alignments` — all genomic alignments for a single query name.  They must
///   all share the same [`query_name`](GenomicAlignment::query_name).  Mixing
///   different names is not an error but will produce incorrect NH/HI values.
/// * `index` — the genome-to-transcriptome index built from your GTF/GFF
///   annotation via [`build_g2t`](crate::g2t::build_g2t).
/// * `config` — evaluation parameters (read type, FASTA availability).
///
/// # Returns
///
/// A `Vec<ProjectedAlignment>` sorted by `(input_index, transcript_id, hi)`.
/// Returns an empty `Vec` when no alignment passes the similarity filter.
pub fn project_group(
    alignments: &[GenomicAlignment],
    index: &G2TTree,
    config: &ProjectionConfig,
) -> Vec<ProjectedAlignment> {
    let evaluator = ReadEvaluator { long_reads: config.long_reads, use_fasta: config.use_fasta };

    let shared_seq: Option<Vec<u8>> = alignments.iter().find_map(|a| a.sequence.clone());

    let mut read_evals: Vec<ReadEval> = Vec::new();

    for (idx, a) in alignments.iter().enumerate() {
        if a.ref_id < 0 { continue; }

        let cigar_ops: Vec<(u32, CigarKind)> = a.cigar.iter()
            .filter_map(|&(len, op)| sam_op_to_kind(op).map(|k| (len, k)))
            .collect();

        let segs = segs_from_ops(a.ref_start, &cigar_ops);
        if segs.is_empty() { continue; }

        let (strand, strand_from_tag) = infer_strand(a);

        let read = ReadAln {
            strand,
            strand_from_tag,
            refid: a.ref_id,
            nh: 0,
            segs,
            cigar_ops,
            sequence: a.sequence.clone(),
            name: a.query_name.clone(),
        };

        let matches: HashMap<Tid, ExonChainMatch> = evaluator
            .evaluate(&read, idx as ReadId, index, shared_seq.as_deref())
            .into_iter()
            .filter(|(_, m)| m.align.cigar.is_some())
            .collect();

        let read_len = if a.read_len > 0 { a.read_len }
            else { a.sequence.as_ref().map_or(0, |s| s.len()) };

        let mut flags = Flags::default();
        if a.is_paired {
            flags.insert(Flags::SEGMENTED);
            if a.mate_is_unmapped { flags.insert(Flags::MATE_UNMAPPED); }
            if a.is_first_in_pair {
                flags.insert(Flags::FIRST_SEGMENT);
            } else {
                flags.insert(Flags::LAST_SEGMENT);
            }
        }
        if a.is_reverse { flags.insert(Flags::REVERSE_COMPLEMENTED); }

        read_evals.push(ReadEval {
            record_idx: idx,
            matches,
            read_len,
            flags,
            alignment_start: Some(a.ref_start as u32),
            mate_alignment_start: a.mate_ref_start.map(|s| s as u32),
            ref_id: Some(a.ref_id),
            mate_ref_id: a.mate_ref_id,
            hit_index: a.hit_index,
        });
    }

    if read_evals.is_empty() {
        return Vec::new();
    }

    let pairs = find_mate_pairs(&read_evals);
    let mut paired = vec![false; read_evals.len()];
    let mut output_groups: Vec<Vec<OutputEntry>> = Vec::new();

    for (i, j) in pairs {
        paired[i] = true;
        paired[j] = true;
        let (r_idx, m_idx) = assign_pair_order(i, j, &read_evals);
        output_groups.extend(build_paired_groups(&read_evals[r_idx], &read_evals[m_idx]));
    }
    for (idx, read) in read_evals.iter().enumerate() {
        if !paired[idx] {
            output_groups.extend(build_unpaired_groups(read));
        }
    }

    let new_nh: u32 = output_groups.iter().map(|g| g.len() as u32).sum();
    let mut entries: Vec<OutputEntry> = output_groups.into_iter().flatten().collect();
    entries.sort_by_key(|e| (e.record_idx, e.tid, e.is_first, e.is_last, e.align.align.hit_index));
    assign_hit_indices(&mut entries);

    let mut out = Vec::with_capacity(entries.len());
    for entry in &entries {
        let is_read1 = entry.mate.is_none() || entry.is_first;
        let nh = if is_read1 { new_nh } else { entry.nh };
        let insert_size = entry.mate.as_ref().map_or(0, |mate| {
            compute_template_length(
                align_pos(&entry.align), mate.pos,
                entry.read_len, mate.read_len,
                entry.is_first, entry.same_transcript,
            )
        });
        out.push(ProjectedAlignment {
            transcript_id: entry.tid,
            transcript_start: entry.align.align.fwpos,
            is_reverse: entry.align.align.is_reverse,
            similarity_score: entry.align.align.similarity_score,
            nh,
            hi: entry.align.align.hit_index as u32,
            is_primary: entry.align.align.primary_alignment,
            same_transcript_as_mate: entry.same_transcript,
            insert_size,
            input_index: entry.record_idx,
        });
    }
    out
}

/// Infer SAM strand character from a [`GenomicAlignment`]'s tag fields.
///
/// Returns `('.'`, false)` when no strand tag is present, signalling the
/// evaluator to check both strands.
fn infer_strand(a: &GenomicAlignment) -> (char, bool) {
    if let Some(s) = a.xs_strand {
        if s == '+' || s == '-' { return (s, true); }
    }
    if let Some(s) = a.ts_strand {
        if s == '+' || s == '-' {
            let flipped = if a.is_reverse { if s == '+' { '-' } else { '+' } } else { s };
            return (flipped, true);
        }
    }
    // No strand tag present — return '.' so the evaluator checks both strands,
    // matching C++ get_strand() behavior when no --fr/--rf flag is set.
    ('.', false)
}
