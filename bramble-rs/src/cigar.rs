//! Pure (htslib-free) projected-CIGAR construction: merge the read's genomic
//! CIGAR with the ideal exon-chain CIGAR into a transcriptome-space SAM CIGAR.
//! Split out of the old `pipeline` module; used by the `bramble-cli` binary's
//! BAM writer. `#[doc(hidden)]` — not part of the stable public API.
#![allow(clippy::all)]
use crate::evaluate::CigarOp;
use noodles::sam::alignment::record::cigar::{op::Kind as CigarKind, Op as SamCigarOp};
use noodles::sam::alignment::record_buf::Cigar as SamCigar;

#[doc(hidden)]
pub fn update_cigar_ops(real_ops: &[SamCigarOp], ideal: &crate::evaluate::Cigar) -> (SamCigar, i32) {
    let (front_hard, front_soft) = front_clip_lengths(real_ops);

    let real_expanded = expand_real_cigar(real_ops);
    let ideal_runs = ideal_runs(ideal);
    let mut ideal_expanded = expand_runs(&ideal_runs);  // mut for pad_ideal_for_leading_clips

    pad_ideal_for_leading_clips(&mut ideal_expanded, front_hard, front_soft);
    let merged = align_and_merge(&real_expanded, &ideal_expanded);

    let mut nm = 0i32;
    let compressed = compress_cigar(merged, &mut nm);
    let cigar = runs_to_sam_cigar(&compressed);

    (cigar, nm)
}

#[allow(dead_code)]
#[doc(hidden)]
pub fn update_cigar_for_test(
    real_ops: Vec<SamCigarOp>,
    ideal: crate::evaluate::Cigar,
) -> (SamCigar, i32) {
    update_cigar_ops(&real_ops, &ideal)
}

fn front_clip_lengths(ops: &[SamCigarOp]) -> (u32, u32) {
    let mut hard = 0u32;
    let mut soft = 0u32;
    let mut idx = 0usize;

    if let Some(op) = ops.first()
        && op.kind() == CigarKind::HardClip
    {
        hard = op.len() as u32;
        idx = 1;
    }

    if let Some(op) = ops.get(idx)
        && op.kind() == CigarKind::SoftClip
    {
        soft = op.len() as u32;
    }

    (hard, soft)
}

fn expand_real_cigar(ops: &[SamCigarOp]) -> Vec<u8> {
    let total_len: usize = ops.iter()
        .filter(|op| op.kind() != CigarKind::Skip)
        .map(|op| op.len())
        .sum();
    let mut expanded = Vec::with_capacity(total_len);
    for op in ops {
        if op.kind() == CigarKind::Skip {
            continue;
        }
        let ch = cigar_kind_to_byte(op.kind());
        expanded.extend(std::iter::repeat_n(ch, op.len()));
    }
    expanded
}

fn cigar_kind_to_byte(kind: CigarKind) -> u8 {
    match kind {
        CigarKind::Match => b'M',
        CigarKind::Insertion => b'I',
        CigarKind::Deletion => b'D',
        CigarKind::Skip => b'N',
        CigarKind::SoftClip => b'S',
        CigarKind::HardClip => b'H',
        CigarKind::Pad => b'P',
        CigarKind::SequenceMatch => b'=',
        CigarKind::SequenceMismatch => b'X',
    }
}

fn ideal_runs(cigar: &crate::evaluate::Cigar) -> Vec<(u32, u8)> {
    let mut runs = Vec::with_capacity(cigar.ops.len());
    for (len, op) in cigar.ops.iter() {
        if *len == 0 {
            continue;
        }
        let ch = match op {
            CigarOp::Match => b'M',
            CigarOp::Ins => b'I',
            CigarOp::Del => b'D',
            CigarOp::RefSkip => b'N',
            CigarOp::SoftClip => b'S',
            CigarOp::HardClip => b'H',
            CigarOp::Pad => b'P',
            CigarOp::Equal => b'=',
            CigarOp::Diff => b'X',
            CigarOp::MatchOverride => b',',
            CigarOp::DelOverride => b'.',
            CigarOp::InsOverride => b'/',
            CigarOp::ClipOverride => b';',
        };
        runs.push((*len, ch));
    }
    runs
}



fn expand_runs(runs: &[(u32, u8)]) -> Vec<u8> {
    let total_len: usize = runs.iter().map(|(len, _)| *len as usize).sum();
    let mut expanded = Vec::with_capacity(total_len);
    for (len, op) in runs {
        expanded.extend(std::iter::repeat_n(*op, *len as usize));
    }
    expanded
}

fn pad_ideal_for_leading_clips(ideal: &mut Vec<u8>, front_hard: u32, front_soft: u32) {
    let total_padding = front_hard.saturating_add(front_soft);
    if total_padding == 0 {
        return;
    }

    let mut padded = Vec::with_capacity(total_padding as usize + ideal.len());

    padded.extend(std::iter::repeat_n(b'_', front_hard as usize));

    // For each soft-clip position, consume a leading override op from ideal
    // (rescue match/del/ins/clip) if one is present, otherwise emit `_` padding.
    // This ensures rescue ops appear BEFORE padding in the merged output,
    // matching C++ merge_cigars which processes override ops first.
    let mut ideal_pos = 0usize;
    for _ in front_hard..front_soft {
        if ideal_pos < ideal.len() && matches!(ideal[ideal_pos], b',' | b'.' | b'/' | b';') {
            padded.push(ideal[ideal_pos]);
            ideal_pos += 1;
        } else {
            padded.push(b'_');
        }
    }

    padded.extend_from_slice(&ideal[ideal_pos..]);
    *ideal = padded;
}

/// Align the expanded real and ideal CIGARs and merge in a single pass,
/// avoiding two intermediate `Vec<u8>` allocations.
fn align_and_merge(real: &[u8], ideal: &[u8]) -> Vec<u8> {
    let max_len = real.len() + ideal.len();
    let mut merged = Vec::with_capacity(max_len);

    let mut real_pos = 0usize;
    let mut ideal_pos = 0usize;
    let max_iterations = max_len * 2;
    let mut iterations = 0usize;

    while real_pos < real.len() || ideal_pos < ideal.len() {
        iterations += 1;
        if iterations > max_iterations {
            break;
        }

        let (r, i);
        if real_pos >= real.len() {
            r = b'_';
            i = ideal[ideal_pos];
            ideal_pos += 1;
        } else if ideal_pos >= ideal.len() {
            r = real[real_pos];
            i = b'_';
            real_pos += 1;
        } else {
            let rb = real[real_pos];
            let ib = ideal[ideal_pos];

            if ib == b'.' {
                r = b'_';
                i = ib;
                ideal_pos += 1;
            } else if rb == b'I' {
                r = rb;
                i = b'_';
                real_pos += 1;
            } else if ib == b'D' {
                r = b'_';
                i = ib;
                ideal_pos += 1;
            } else {
                r = rb;
                i = ib;
                real_pos += 1;
                ideal_pos += 1;
            }
        }

        let merge = merge_ops(r, i);
        // if a _ is returned, nothing should be added
        if merge != b'_' {
            merged.push(merge);   
        } 
    }

    merged
}

#[inline(always)]
fn merge_ops(real_op: u8, ideal_op: u8) -> u8 {
    // Override ops from clip rescue — handle first as they dominate in FASTA mode.
    match ideal_op {
        b';' => return if real_op == b'D' { b'_' } else { b'S' },
        b',' => return if real_op == b'D' { b'D' } else if real_op == b'I' { b'I' } else { b'M' },
        b'/' => return if real_op == b'D' || real_op == b'I' { b'_' } else { b'I' },
        b'.' => return if real_op == b'D' || real_op == b'I' { b'_' } else { b'D' },
        b'*' => return real_op,
        b'_' => return real_op,
        _ => {}
    }

    // real_op == I with non-override ideal: special-case _
    if real_op == b'I' && ideal_op == b'_' {
        return b'I';
    }

    if real_op == b'H' {
        return b'H';
    }

    if real_op == b'P' {
        return ideal_op;
    }

    // Standard CIGAR ops
    match (real_op, ideal_op) {
        (b'D', b'S') => b'_',
        (b'I', b'S') => b'S',
        (b'D', b'I') => b'_',
        _ if ideal_op == b'S' || ideal_op == b'D' || ideal_op == b'I' => ideal_op,
        _ if real_op == b'S' || real_op == b'D' || real_op == b'I' => real_op,
        _ if ideal_op == b'M' || ideal_op == b'=' || ideal_op == b'X' => b'M',
        _ if real_op == b'M' || real_op == b'=' || real_op == b'X' => b'M',
        (b'_', _) => ideal_op,
        (_, b'_') => real_op,
        _ => if real_op != b'*' { real_op } else { ideal_op },
    }
}

fn compress_cigar(cleaned: Vec<u8>, nm: &mut i32) -> Vec<(u32, u8)> {
    let mut runs = Vec::new();
    let mut i = 0usize;
    while i < cleaned.len() {
        if cleaned[i] == b'_' {
            i += 1;
            continue;
        }

        let op_byte = cleaned[i];
        let mut run_len = 0u32;
        while i < cleaned.len() && (cleaned[i] == op_byte || cleaned[i] == b'_') {
            if cleaned[i] == op_byte {
                run_len += 1;
            }
            i += 1;
        }

        if run_len == 0 {
            continue;
        }

        let out_byte = match op_byte {
            b'M' | b'=' | b'X' => b'M',
            b'I' => {
                *nm += run_len as i32;
                b'I'
            }
            b'D' => {
                *nm += run_len as i32;
                b'D'
            }
            b'S' => b'S',
            b'H' => b'H',
            b'N' => b'N',
            b'P' => b'P',
            _ => continue,
        };

        runs.push((run_len, out_byte));
    }

    runs
}

fn runs_to_sam_cigar(runs: &[(u32, u8)]) -> SamCigar {
    let ops: Vec<SamCigarOp> = runs
        .iter()
        .filter_map(|(len, op)| {
            if *len == 0 {
                return None;
            }
            let kind = match op {
                b'M' => CigarKind::Match,
                b'I' => CigarKind::Insertion,
                b'D' => CigarKind::Deletion,
                b'N' => CigarKind::Skip,
                b'S' => CigarKind::SoftClip,
                b'H' => CigarKind::HardClip,
                b'P' => CigarKind::Pad,
                b'=' => CigarKind::SequenceMatch,
                b'X' => CigarKind::SequenceMismatch,
                _ => return None,
            };
            Some(SamCigarOp::new(kind, *len as usize))
        })
        .collect();

    ops.into_iter().collect()
}
