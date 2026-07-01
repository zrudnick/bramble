//! Library-internal helpers shared by the projection API ([`crate::api`]) and the
//! `bramble-cli` binary: read/output grouping, mate pairing, hit-index assignment,
//! and CIGAR-op conversion. Split out of the old `pipeline` module so the library
//! no longer depends on rust-htslib (BAM I/O lives in `bramble-cli`). Items are
//! `pub` so the binary crate can use them, but `#[doc(hidden)]` (not stable API).
#![allow(clippy::all)]
use crate::alignment;
use crate::evaluate::ExonChainMatch;
use crate::types::{HashMap, HashMapExt, HashSet, RefId, Tid};
use noodles::sam::alignment::record::cigar::op::Kind as CigarKind;

#[derive(Debug)]
#[doc(hidden)]
pub struct ReadEval {
    pub record_idx: usize,
    pub matches: HashMap<Tid, ExonChainMatch>,
    pub read_len: usize,
    pub flags: u16,
    pub alignment_start: Option<u32>,
    pub mate_alignment_start: Option<u32>,
    pub ref_id: Option<RefId>,
    pub mate_ref_id: Option<RefId>,
    pub hit_index: i32,
}

// BAM flag constants
#[doc(hidden)]
pub const FLAG_PAIRED: u16       = 0x1;
#[doc(hidden)]
pub const FLAG_PROPER_PAIR: u16  = 0x2;
#[doc(hidden)]
pub const FLAG_UNMAPPED: u16     = 0x4;
#[doc(hidden)]
pub const FLAG_MATE_UNMAPPED: u16 = 0x8;
#[doc(hidden)]
pub const FLAG_REVERSE: u16      = 0x10;
#[doc(hidden)]
pub const FLAG_MATE_REVERSE: u16 = 0x20;
#[doc(hidden)]
pub const FLAG_READ1: u16        = 0x40;
#[doc(hidden)]
pub const FLAG_READ2: u16        = 0x80;
#[doc(hidden)]
pub const FLAG_SECONDARY: u16    = 0x100;

#[derive(Debug, Clone)]
#[doc(hidden)]
pub struct MateInfo {
    pub tid: Tid,
    pub pos: u32,
    pub is_reverse: bool,
}

#[derive(Debug, Clone)]
#[doc(hidden)]
pub struct OutputEntry {
    pub record_idx: usize,
    pub tid: Tid,
    pub align: ExonChainMatch,
    pub read_len: usize,
    pub mate: Option<MateInfo>,
    pub is_first: bool,
    pub is_last: bool,
    pub same_transcript: bool,
    pub nh: u32,
}

#[doc(hidden)]
pub fn assign_hit_indices(entries: &mut [OutputEntry]) {
    let mut by_read: HashMap<usize, Vec<usize>> = HashMap::new();
    for (idx, entry) in entries.iter().enumerate() {
        by_read.entry(entry.record_idx).or_default().push(idx);
    }

    for (_record_idx, mut indices) in by_read {
        indices.sort_by(|&a, &b| {
            entries[a]
                .tid
                .cmp(&entries[b].tid)
                .then_with(|| entries[a].is_first.cmp(&entries[b].is_first))
                .then_with(|| entries[a].is_last.cmp(&entries[b].is_last))
        });

        let mut best_score = f64::NEG_INFINITY;
        let mut best_idx: Option<usize> = None;

        for (hit_index, entry_idx) in indices.iter().copied().enumerate() {
            let entry = &mut entries[entry_idx];
            entry.align.align.hit_index = (hit_index + 1) as i32;
            let score = entry.align.align.similarity_score;
            if score > best_score + 1e-12 {
                best_score = score;
                best_idx = Some(entry_idx);
            } else if (score - best_score).abs() <= 1e-12
                && let Some(current) = best_idx
                && entry.tid < entries[current].tid
            {
                best_idx = Some(entry_idx);
            }
        }

        let best_idx = best_idx.unwrap_or(indices[0]);
        for entry_idx in indices.iter().copied() {
            entries[entry_idx].align.align.primary_alignment = entry_idx == best_idx;
        }
    }

    // C++ marks BOTH mates in the best pair as primary. Propagate primary
    // status: if a paired entry is primary, find its mate (same tid, opposite
    // is_first/is_last) and mark it primary too.
    let primary_tids: Vec<(Tid, bool)> = entries
        .iter()
        .filter(|e| e.align.align.primary_alignment && e.mate.is_some())
        .map(|e| (e.tid, e.is_first))
        .collect();
    for (tid, is_first) in primary_tids {
        for entry in entries.iter_mut() {
            if entry.mate.is_some() && entry.tid == tid && entry.is_first != is_first {
                entry.align.align.primary_alignment = true;
            }
        }
    }
}


#[doc(hidden)]
pub fn find_mate_pairs(reads: &[ReadEval]) -> Vec<(usize, usize)> {
    // Pair mates that MUTUALLY reference each other's alignment start (same name
    // group, same reference, read.mate_start == mate.read_start and vice versa),
    // independent of the order reads appear in the group.
    //
    // The previous single-pass logic (mirroring C++ `hashread[key]=id`) only
    // *inserted* the LEFT mate (mate_start > read_start) and only *looked up* from
    // the RIGHT mate, so it paired a fragment ONLY when the LEFT mate was processed
    // before the RIGHT mate. When the RIGHT (higher-coordinate) mate came first —
    // e.g. a name-collated BAM listing read1 first while read1 is the reverse mate
    // — the lookup hit an empty map and the LEFT mate's later insert was never
    // consumed, so BOTH mates fell through to `build_unpaired_groups` and were
    // emitted as orphans on every transcript (losing the proper-pair fragment
    // length signal). C++ bramble has the same order bug. Indexing all eligible
    // reads first, then matching mutually, removes the order dependence.
    let eligible = |r: &ReadEval| -> Option<(u32, u32)> {
        if (r.flags & FLAG_PAIRED) == 0 || (r.flags & FLAG_MATE_UNMAPPED) != 0 {
            return None;
        }
        match (
            r.alignment_start,
            r.mate_alignment_start,
            r.ref_id,
            r.mate_ref_id,
        ) {
            (Some(rs), Some(ms), Some(rid), Some(mid)) if rid == mid => Some((rs, ms)),
            _ => None,
        }
    };

    // Index eligible reads by their own alignment start (insertion order kept for
    // deterministic candidate selection among same-position multimappers).
    let mut by_start: HashMap<u32, Vec<usize>> = HashMap::new();
    for (idx, read) in reads.iter().enumerate() {
        if let Some((rs, _)) = eligible(read) {
            by_start.entry(rs).or_default().push(idx);
        }
    }

    let mut paired = vec![false; reads.len()];
    let mut pairs = Vec::new();
    for (idx, read) in reads.iter().enumerate() {
        if paired[idx] {
            continue;
        }
        let Some((rs, ms)) = eligible(read) else {
            continue;
        };
        // Our mate sits at `ms`; find an unpaired read there that points back at us.
        if let Some(cands) = by_start.get(&ms) {
            for &m in cands {
                if m == idx || paired[m] {
                    continue;
                }
                if let Some((_, m_ms)) = eligible(&reads[m])
                    && m_ms == rs
                {
                    pairs.push((idx, m));
                    paired[idx] = true;
                    paired[m] = true;
                    break;
                }
            }
        }
    }

    pairs
}

#[doc(hidden)]
pub fn assign_pair_order(i: usize, j: usize, reads: &[ReadEval]) -> (usize, usize) {
    let flags_i = reads[i].flags;
    let flags_j = reads[j].flags;

    let i_first = (flags_i & FLAG_READ1) != 0;
    let i_last = (flags_i & FLAG_READ2) != 0;
    let j_first = (flags_j & FLAG_READ1) != 0;
    let j_last = (flags_j & FLAG_READ2) != 0;

    if i_first && j_last {
        return (i, j);
    }
    if j_first && i_last {
        return (j, i);
    }
    if i_first && !j_first {
        return (i, j);
    }
    if j_first && !i_first {
        return (j, i);
    }

    if i <= j {
        (i, j)
    } else {
        (j, i)
    }
}

#[doc(hidden)]
pub fn build_paired_groups(read: &ReadEval, mate: &ReadEval) -> Vec<Vec<OutputEntry>> {
    let read_transcripts: HashSet<Tid> = read.matches.keys().copied().collect();
    let mate_transcripts: HashSet<Tid> = mate.matches.keys().copied().collect();
    let read_nh = read.matches.len() as u32;
    let mate_nh = mate.matches.len() as u32;

    if read_transcripts.is_empty() || mate_transcripts.is_empty() {
        return Vec::new();
    }

    let common: HashSet<Tid> = read_transcripts
        .intersection(&mate_transcripts)
        .copied()
        .collect();

    let mut groups = Vec::new();

    if !common.is_empty() {
        for tid in common {
            let Some(read_match) = read.matches.get(&tid) else { continue; };
            let Some(mate_match) = mate.matches.get(&tid) else { continue; };

            let mate_info_for_read = MateInfo {
                tid,
                pos: align_pos(mate_match),
                is_reverse: mate_match.align.is_reverse,
            };
            let mate_info_for_mate = MateInfo {
                tid,
                pos: align_pos(read_match),
                is_reverse: read_match.align.is_reverse,
            };

            groups.push(vec![
                OutputEntry {
                    record_idx: read.record_idx,
                    tid,
                    align: read_match.clone(),
                    read_len: read.read_len,
                    mate: Some(mate_info_for_read),
                    is_first: true,
                    is_last: false,
                    same_transcript: true,
                    nh: read_nh,
                },
                OutputEntry {
                    record_idx: mate.record_idx,
                    tid,
                    align: mate_match.clone(),
                    read_len: mate.read_len,
                    mate: Some(mate_info_for_mate),
                    is_first: false,
                    is_last: true,
                    same_transcript: true,
                    nh: mate_nh,
                },
            ]);
        }

        return groups;
    }

    if read_transcripts.len() == 1 && mate_transcripts.len() == 1 {
        let r_tid = *read_transcripts.iter().next().unwrap();
        let m_tid = *mate_transcripts.iter().next().unwrap();

        let Some(read_match) = read.matches.get(&r_tid) else { return groups; };
        let Some(mate_match) = mate.matches.get(&m_tid) else { return groups; };

        let mate_info_for_read = MateInfo {
            tid: m_tid,
            pos: align_pos(mate_match),
            is_reverse: mate_match.align.is_reverse,
        };
        let mate_info_for_mate = MateInfo {
            tid: r_tid,
            pos: align_pos(read_match),
            is_reverse: read_match.align.is_reverse,
        };

        groups.push(vec![
            OutputEntry {
                record_idx: read.record_idx,
                tid: r_tid,
                align: read_match.clone(),
                read_len: read.read_len,
                mate: Some(mate_info_for_read),
                is_first: true,
                is_last: false,
                same_transcript: false,
                nh: read_nh,
            },
            OutputEntry {
                record_idx: mate.record_idx,
                tid: m_tid,
                align: mate_match.clone(),
                read_len: mate.read_len,
                mate: Some(mate_info_for_mate),
                is_first: false,
                is_last: true,
                same_transcript: false,
                nh: mate_nh,
            },
        ]);
    }

    groups
}

#[doc(hidden)]
pub fn build_unpaired_groups(read: &ReadEval) -> Vec<Vec<OutputEntry>> {
    let mut groups = Vec::new();
    let read_nh = read.matches.len() as u32;

    for (tid, align) in read.matches.iter() {
        groups.push(vec![OutputEntry {
            record_idx: read.record_idx,
            tid: *tid,
            align: align.clone(),
            read_len: read.read_len,
            mate: None,
            is_first: false,
            is_last: false,
            same_transcript: false,
            nh: read_nh,
        }]);
    }

    groups
}

#[doc(hidden)]
pub fn align_pos(match_info: &ExonChainMatch) -> u32 {
    if match_info.align.strand == '+' {
        match_info.align.fwpos
    } else {
        match_info.align.rcpos
    }
}

#[doc(hidden)]
pub fn compute_template_length(
    my_pos: u32,
    mate_pos: u32,
    read_len: usize,
    same_transcript: bool,
) -> i32 {
    if !same_transcript {
        return 0;
    }

    let my_pos = my_pos as i32;
    let mate_pos = mate_pos as i32;
    let l_qseq = read_len as i32;

    // C++ SAM spec: positive for leftmost read, negative for rightmost
    if my_pos <= mate_pos {
        (mate_pos + l_qseq) - my_pos
    } else {
        -((my_pos + l_qseq) - mate_pos)
    }
}


#[doc(hidden)]
pub fn sam_op_to_kind(op: u8) -> Option<CigarKind> {
    match op {
        0 => Some(CigarKind::Match),
        1 => Some(CigarKind::Insertion),
        2 => Some(CigarKind::Deletion),
        3 => Some(CigarKind::Skip),
        4 => Some(CigarKind::SoftClip),
        5 => Some(CigarKind::HardClip),
        6 => Some(CigarKind::Pad),
        7 => Some(CigarKind::SequenceMatch),
        8 => Some(CigarKind::SequenceMismatch),
        _ => None,
    }
}

/// Compute exon segments from a raw CIGAR op list and a 1-based reference start.
#[doc(hidden)]
pub fn segs_from_ops(ref_start: i64, ops: &[(u32, CigarKind)]) -> Vec<alignment::Segment> {
    let mut ref_pos = ref_start as u32;
    let mut exon_start = ref_pos;
    let mut segs = Vec::new();
    for &(len, kind) in ops {
        match kind {
            CigarKind::Match
            | CigarKind::SequenceMatch
            | CigarKind::SequenceMismatch
            | CigarKind::Deletion => {
                ref_pos = ref_pos.saturating_add(len);
            }
            CigarKind::Skip => {
                if ref_pos > exon_start {
                    segs.push(alignment::Segment { start: exon_start, end: ref_pos });
                }
                ref_pos = ref_pos.saturating_add(len);
                exon_start = ref_pos;
            }
            _ => {}
        }
    }
    if ref_pos > exon_start {
        segs.push(alignment::Segment { start: exon_start, end: ref_pos });
    }
    segs
}


