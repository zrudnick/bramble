//! Pure (htslib-free) projected-CIGAR construction: merge the read's genomic
//! CIGAR with the ideal exon-chain CIGAR into a transcriptome-space SAM CIGAR.
//! Split out of the old `pipeline` module; used by the `bramble-cli` binary's
//! BAM writer. `#[doc(hidden)]` — not part of the stable public API.
#![allow(clippy::all)]
use crate::evaluate::CigarOp;
use noodles::sam::alignment::record::cigar::{Op as SamCigarOp, op::Kind as CigarKind};
use noodles::sam::alignment::record_buf::Cigar as SamCigar;

#[doc(hidden)]
pub fn update_cigar_ops(
    real_ops: &[SamCigarOp],
    ideal: &crate::evaluate::Cigar,
) -> (SamCigar, i32) {
    let (front_hard, front_soft) = front_clip_lengths(real_ops);

    let real_expanded = expand_real_cigar(real_ops);
    let ideal_runs = ideal_runs(ideal);
    let mut ideal_expanded = expand_runs(&ideal_runs); // mut for pad_ideal_for_leading_clips

    pad_ideal_for_leading_clips(&mut ideal_expanded, front_hard, front_soft);
    let merged = align_and_merge(&real_expanded, &ideal_expanded);

    let mut nm = 0i32;
    let mut compressed = compress_cigar(merged, &mut nm);
    absorb_insertions_between_clips(&mut compressed, &mut nm);
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
    let total_len: usize = ops
        .iter()
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
    for _ in 0..front_soft {
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
        b',' => {
            return if real_op == b'D' {
                b'D'
            } else if real_op == b'I' {
                b'I'
            } else {
                b'M'
            };
        }
        b'/' => {
            return if real_op == b'D' || real_op == b'I' {
                b'_'
            } else {
                b'I'
            };
        }
        b'.' => {
            return if real_op == b'D' || real_op == b'I' {
                b'_'
            } else {
                b'D'
            };
        }
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
        _ => {
            if real_op != b'*' {
                real_op
            } else {
                ideal_op
            }
        }
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

/// An insertion between two clips (e.g. `9S 1I 17S`, from a read insertion that
/// falls in the clipped overhang past the transcript end) is part of the clip:
/// turn it into the clip op and merge the neighbours, giving `27S`. A clip in
/// the middle of a CIGAR is invalid SAM. Mirrors the C++ pass in merge_cigars
/// (bam.cpp). The insertion's bases are no longer edits, so they leave `nm`.
fn absorb_insertions_between_clips(runs: &mut Vec<(u32, u8)>, nm: &mut i32) {
    let is_clip = |op: u8| op == b'S' || op == b'H';
    for i in 1..runs.len().saturating_sub(1) {
        if runs[i].1 == b'I' && is_clip(runs[i - 1].1) && is_clip(runs[i + 1].1) {
            *nm -= runs[i].0 as i32;
            runs[i].1 = runs[i - 1].1;
        }
    }
    // merge neighbouring runs that now share an op
    let mut merged: Vec<(u32, u8)> = Vec::with_capacity(runs.len());
    for &(len, op) in runs.iter() {
        match merged.last_mut() {
            Some(last) if last.1 == op => last.0 += len,
            _ => merged.push((len, op)),
        }
    }
    *runs = merged;
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

/// Count the MD mismatches that fall in `keep` (sorted, disjoint query
/// intervals [start, end), genome orientation), walking the MD tag and the
/// genome CIGAR together without allocating (this runs once per output record).
/// `None` if MD is malformed or does not fit the CIGAR. Mirrors C++
/// count_md_mismatches (bam.cpp).
fn count_md_mismatches(genome_ops: &[SamCigarOp], md: &[u8], keep: &[(u32, u32)]) -> Option<i64> {
    let is_aligned =
        |k: CigarKind| matches!(k, CigarKind::Match | CigarKind::SequenceMatch | CigarKind::SequenceMismatch);
    let mut ci = 0usize; // current genome CIGAR op
    let mut off = 0u32; // aligned bases already used in op ci
    let mut q = 0u32; // query position at the start of op ci
    // move to the next unused reference-aligned base; false if none left
    let next_aligned = |ci: &mut usize, off: &mut u32, q: &mut u32| -> bool {
        while *ci < genome_ops.len() {
            let op = &genome_ops[*ci];
            let len = op.len() as u32;
            if is_aligned(op.kind()) && *off < len {
                return true;
            }
            if op.kind().consumes_read() {
                *q += len;
            }
            *ci += 1;
            *off = 0;
        }
        false
    };

    let mut count = 0i64;
    let mut kp = 0usize;
    let mut i = 0usize;
    while i < md.len() {
        let c = md[i];
        if c.is_ascii_digit() {
            let mut n = 0u64;
            while i < md.len() && md[i].is_ascii_digit() {
                n = n * 10 + (md[i] - b'0') as u64;
                i += 1;
            }
            while n > 0 {
                // matching bases
                if !next_aligned(&mut ci, &mut off, &mut q) {
                    return None;
                }
                let avail = genome_ops[ci].len() as u32 - off;
                let take = (n.min(avail as u64)) as u32;
                off += take;
                n -= take as u64;
            }
        } else if c == b'^' {
            // deleted reference bases: no query position
            i += 1;
            while i < md.len() && md[i].is_ascii_alphabetic() {
                i += 1;
            }
        } else if c.is_ascii_alphabetic() {
            if !next_aligned(&mut ci, &mut off, &mut q) {
                return None;
            }
            let pos = q + off;
            off += 1;
            while kp < keep.len() && keep[kp].1 <= pos {
                kp += 1;
            }
            if kp < keep.len() && keep[kp].0 <= pos {
                count += 1;
            }
            i += 1;
        } else {
            return None;
        }
    }
    // MD must cover every aligned base
    if next_aligned(&mut ci, &mut off, &mut q) { None } else { Some(count) }
}

thread_local! {
    // aligned query intervals of the projected CIGAR, reused across records
    static KEEP: std::cell::RefCell<Vec<(u32, u32)>> = const { std::cell::RefCell::new(Vec::new()) };
}

/// NM for the projected alignment: mismatches of bases still aligned outside
/// clip rescue (from MD; else genome NM minus genome indels), plus mismatches in
/// rescued bases, plus every inserted/deleted base of the projected CIGAR.
/// `projected_ops` must be in genome orientation (before reversing for '-').
/// `None` when the genome record has neither NM nor MD. Mirrors C++
/// set_projected_nm (bam.cpp).
#[doc(hidden)]
pub fn projected_nm(
    genome_ops: &[SamCigarOp],
    projected_ops: &[SamCigarOp],
    genome_nm: Option<i64>,
    md: Option<&[u8]>,
    rescue_left: u32,
    rescue_right: u32,
    rescue_mm: u32,
) -> Option<i32> {
    if genome_nm.is_none() && md.is_none() {
        return None;
    }
    let qlen: u32 = projected_ops
        .iter()
        .filter(|op| op.kind().consumes_read())
        .map(|op| op.len() as u32)
        .sum();
    let lo = rescue_left.min(qlen);
    let hi = qlen - rescue_right.min(qlen);

    KEEP.with(|keep| {
        // aligned (M/=/X) query intervals of the projected CIGAR, outside the
        // ranges realigned by clip rescue
        let mut keep = keep.borrow_mut();
        keep.clear();
        let mut new_indels = 0i64;
        let mut q = 0u32;
        for op in projected_ops {
            let len = op.len() as u32;
            match op.kind() {
                CigarKind::Match | CigarKind::SequenceMatch | CigarKind::SequenceMismatch => {
                    let (s, e) = (q.max(lo), (q + len).min(hi));
                    if s < e {
                        keep.push((s, e));
                    }
                }
                CigarKind::Insertion | CigarKind::Deletion => new_indels += len as i64,
                _ => {}
            }
            if op.kind().consumes_read() {
                q += len;
            }
        }

        let mismatches = match md.and_then(|md| count_md_mismatches(genome_ops, md, &keep)) {
            Some(m) => m,
            None => {
                // no usable MD
                let genome_nm = genome_nm?;
                let genome_indels: i64 = genome_ops
                    .iter()
                    .filter(|op| matches!(op.kind(), CigarKind::Insertion | CigarKind::Deletion))
                    .map(|op| op.len() as i64)
                    .sum();
                (genome_nm - genome_indels).max(0)
            }
        };
        Some((mismatches + rescue_mm as i64 + new_indels) as i32)
    })
}
