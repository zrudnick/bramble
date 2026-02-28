//! Public library API for projecting spliced genomic alignments into transcriptome space.
//!
//! # Example
//!
//! ```no_run
//! use bramble_rs::{GenomicAlignment, ProjectionConfig, project_group};
//! use bramble_rs::g2t::G2TTree;
//!
//! // Build the index from a GTF/GFF annotation:
//! // let transcripts = bramble_rs::annotation::load_transcripts(path)?;
//! // let index = bramble_rs::g2t::build_g2t(&transcripts, &refname_to_id, None)?;
//! //
//! // let config = ProjectionConfig { long_reads: false, use_fasta: false };
//! // let alns: Vec<GenomicAlignment> = /* from minimap2-rs or noodles */;
//! // let projected = project_group(&alns, &index, &config);
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
/// Construct from noodles `RecordBuf` or from minimap2-rs `Mapping`.
#[derive(Debug, Clone)]
pub struct GenomicAlignment {
    pub query_name: String,
    /// 0-based reference sequence index (matches the index used to build `G2TTree`).
    pub ref_id: i32,
    /// 1-based alignment start position on the reference.
    pub ref_start: i64,
    /// True if the read maps to the reverse strand.
    pub is_reverse: bool,
    /// CIGAR as `(length, SAM op code)` pairs.
    /// SAM op codes: 0=M, 1=I, 2=D, 3=N(intron), 4=S, 5=H, 6=P, 7==, 8=X.
    pub cigar: Vec<(u32, u8)>,
    /// Query sequence bytes (ASCII bases). Required for soft-clip rescue.
    pub sequence: Option<Vec<u8>>,
    /// True if the read comes from a paired-end library.
    pub is_paired: bool,
    /// True if this is read1 of a pair (FIRST_SEGMENT flag in SAM).
    pub is_first_in_pair: bool,
    /// Library strand inferred from the XS tag (short reads).
    pub xs_strand: Option<char>,
    /// Library strand inferred from the ts tag (long reads).
    pub ts_strand: Option<char>,
    /// HI tag value from the input (0 if absent).
    pub hit_index: i32,
    /// Reference sequence index of the mate (None if unpaired or mate unmapped).
    pub mate_ref_id: Option<i32>,
    /// 1-based alignment start of the mate (None if unpaired or mate unmapped).
    pub mate_ref_start: Option<i64>,
    /// True if the mate is unmapped.
    pub mate_is_unmapped: bool,
    /// Query read length in bases (used for NH when sequence is absent).
    pub read_len: usize,
}

/// The result of projecting one alignment into transcriptome space.
#[derive(Debug, Clone)]
pub struct ProjectedAlignment {
    /// Index of the target transcript in the `G2TTree`.
    pub transcript_id: u32,
    /// 1-based start position on the transcript.
    pub transcript_start: u32,
    /// True if the alignment is on the reverse strand of the transcript.
    pub is_reverse: bool,
    /// Alignment similarity score (higher is better).
    pub similarity_score: f64,
    /// Total hits for this read (NH tag value).
    pub nh: u32,
    /// 1-based index of this hit (HI tag value).
    pub hi: u32,
    /// True if this is the primary (best) alignment.
    pub is_primary: bool,
    /// True if the mate maps to the same transcript.
    pub same_transcript_as_mate: bool,
    /// Observed template length (TLEN); 0 when mates are on different transcripts.
    pub insert_size: i32,
    /// The index of the input `GenomicAlignment` this result was derived from.
    pub input_index: usize,
}

/// Configuration passed to `project_group()`.
#[derive(Debug, Clone)]
pub struct ProjectionConfig {
    /// Use long-read evaluation parameters (more permissive clip/gap tolerances).
    pub long_reads: bool,
    /// Whether a FASTA sequence is available for clip-rescue alignment.
    pub use_fasta: bool,
}

/// Project a group of alignments for the same query name into transcriptome coordinates.
///
/// All alignments in the slice must share the same `query_name`. The returned
/// `ProjectedAlignment` slice is sorted by `(input_index, transcript_id, hi)`.
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

/// Infer SAM strand character from a `GenomicAlignment`'s tag fields.
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
