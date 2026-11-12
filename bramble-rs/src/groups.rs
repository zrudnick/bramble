//! Library-internal helpers shared by the projection API ([`crate::api`]) and the
//! `bramble-cli` binary: read/output grouping, mate pairing, hit-index assignment,
//! and CIGAR-op conversion. Split out of the old `pipeline` module so the library
//! no longer depends on rust-htslib (BAM I/O lives in `bramble-cli`). Items are
//! `pub` so the binary crate can use them, but `#[doc(hidden)]` (not stable API).
#![allow(clippy::all)]
use crate::alignment;
use crate::evaluate::ExonChainMatch;
use crate::g2t::G2TTree;
use crate::types::{HashMap, HashMapExt, HashSet, HashSetExt, RefId, Tid};
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
pub const FLAG_PAIRED: u16 = 0x1;
#[doc(hidden)]
pub const FLAG_PROPER_PAIR: u16 = 0x2;
#[doc(hidden)]
pub const FLAG_UNMAPPED: u16 = 0x4;
#[doc(hidden)]
pub const FLAG_MATE_UNMAPPED: u16 = 0x8;
#[doc(hidden)]
pub const FLAG_REVERSE: u16 = 0x10;
#[doc(hidden)]
pub const FLAG_MATE_REVERSE: u16 = 0x20;
#[doc(hidden)]
pub const FLAG_READ1: u16 = 0x40;
#[doc(hidden)]
pub const FLAG_READ2: u16 = 0x80;
#[doc(hidden)]
pub const FLAG_SECONDARY: u16 = 0x100;

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

/// Collapse duplicate output units of one read. A record is identified by its
/// genomic alignment (`genome_key` of its record index: read end, strand,
/// reference, position, CIGAR) and its transcript placement; the same genomic
/// alignment can be reported once per mate (e.g. three identical read1
/// alignments with three different mates), giving the same transcript record
/// several times. Pairs are kept first, then single records not already kept,
/// each in input order, so the earliest copy survives. Mirrors C++
/// collapse_duplicates (core.cpp). Call before NH, HI and the primary.
#[doc(hidden)]
pub fn collapse_duplicate_units<K: std::hash::Hash + Eq + Clone>(
    groups: &mut Vec<Vec<OutputEntry>>,
    genome_key: impl Fn(usize) -> K,
) {
    let rec_key = |e: &OutputEntry| (genome_key(e.record_idx), e.tid, align_pos(&e.align), e.align.align.strand);
    let mut units: Vec<Vec<OutputEntry>> = std::mem::take(groups).into_iter().filter(|g| !g.is_empty()).collect();
    // input order: pairs by (read, mate) record index, then singles by record index
    units.sort_by_key(|g| (g.len() != 2, g[0].record_idx, g.get(1).map_or(0, |m| m.record_idx), g[0].tid));

    let mut kept_records = HashSet::new();
    let mut kept_pairs = HashSet::new();
    for g in units {
        if g.len() == 2 {
            let (a, b) = (rec_key(&g[0]), rec_key(&g[1]));
            if kept_pairs.insert((a.clone(), b.clone())) {
                kept_records.insert(a);
                kept_records.insert(b);
                groups.push(g);
            }
        } else if g.iter().all(|e| !kept_records.contains(&rec_key(e))) {
            for e in &g {
                kept_records.insert(rec_key(e));
            }
            groups.push(g);
        }
    }
}

/// Rank every projected alignment of one read together, as C++ does when it
/// flushes a read-name group (core.cpp). Each output group is one unit: a mate
/// pair (read, then mate) or a single unpaired record. HI runs 1..NH across the
/// whole read, the unit with the highest score (max over its records) is
/// primary (both mates for a pair), and ties are broken by a hash of the read
/// name so the choice does not depend on transcript id or processing order.
/// Tied units are sorted by transcript name and position first, matching C++.
#[doc(hidden)]
pub fn assign_primary_and_hit_indices(
    groups: &mut [Vec<OutputEntry>],
    name: &[u8],
    g2t: &G2TTree,
    read_end: impl Fn(usize) -> u8,
) {
    // deterministic unit order, independent of HashMap iteration order
    groups.sort_by_key(|g| g.first().map(|e| (e.record_idx, e.tid)));

    let mut hit_index: i32 = 0;
    let mut best_score = f64::NEG_INFINITY;
    let mut tied: Vec<usize> = Vec::new();

    for (unit, group) in groups.iter_mut().enumerate() {
        if group.is_empty() {
            continue;
        }
        let mut unit_score = f64::NEG_INFINITY;
        for entry in group.iter_mut() {
            hit_index += 1;
            entry.align.align.hit_index = hit_index;
            entry.align.align.primary_alignment = false;
            unit_score = unit_score.max(entry.align.align.similarity_score);
        }
        if unit_score > best_score {
            best_score = unit_score;
            tied.clear();
            tied.push(unit);
        } else if unit_score == best_score {
            tied.push(unit);
        }
    }

    if tied.is_empty() {
        return;
    }
    // sort ties by transcript names and positions, not tids: tid numbering
    // differs between the C++ and Rust versions (C++ core.cpp tie_key)
    let tie_key = |g: &Vec<OutputEntry>| {
        let tname = |e: &OutputEntry| g2t.transcript_name(e.tid).unwrap_or("").to_string();
        // then input order of the genomic alignments (same in C++), for ties on
        // the same transcript and position
        let r = &g[0];
        match g.get(1) {
            None => (tname(r), align_pos(&r.align), false, String::new(), 0, r.record_idx, 0),
            Some(m) => (
                tname(r),
                align_pos(&r.align),
                true,
                tname(m),
                align_pos(&m.align),
                r.record_idx,
                m.record_idx,
            ),
        }
    };
    tied.sort_by_key(|&u| tie_key(&groups[u]));
    let best = tied[crate::evaluate::deterministic_index(name, tied.len())];
    for entry in groups[best].iter_mut() {
        entry.align.align.primary_alignment = true;
    }

    // Every read end gets one primary. If the best unit is a single record, the
    // other end's best record (any unit) becomes primary on its own; its partner,
    // if any, stays secondary. Ties: sorted by transcript name, position and input
    // order, picked by the read-name hash (as C++ core.cpp).
    let covered: Vec<u8> = groups[best].iter().map(|e| read_end(e.record_idx)).collect();
    let mut ends: Vec<u8> = groups.iter().flatten().map(|e| read_end(e.record_idx)).collect();
    ends.sort_unstable();
    ends.dedup();
    for end in ends {
        if end == 0 || covered.contains(&end) {
            continue;
        }
        let mut best_score = f64::NEG_INFINITY;
        let mut cands: Vec<(usize, usize)> = Vec::new();
        for (u, g) in groups.iter().enumerate() {
            for (k, e) in g.iter().enumerate() {
                if read_end(e.record_idx) != end {
                    continue;
                }
                let score = e.align.align.similarity_score;
                if score > best_score {
                    best_score = score;
                    cands.clear();
                }
                if score == best_score {
                    cands.push((u, k));
                }
            }
        }
        if cands.is_empty() {
            continue;
        }
        let rec_key = |&(u, k): &(usize, usize)| {
            let e = &groups[u][k];
            let tname = g2t.transcript_name(e.tid).unwrap_or("").to_string();
            (tname, align_pos(&e.align), e.record_idx)
        };
        cands.sort_by_key(rec_key);
        let (u, k) = cands[crate::evaluate::deterministic_index(name, cands.len())];
        groups[u][k].align.align.primary_alignment = true;
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
                // opposite read ends only: read1 pairs with read2
                let end = |r: &ReadEval| r.flags & (FLAG_READ1 | FLAG_READ2);
                if let Some((_, m_ms)) = eligible(&reads[m])
                    && m_ms == rs
                    && end(&reads[m]) != end(read)
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

    if i <= j { (i, j) } else { (j, i) }
}

/// Orientation of a record in the output: its genome orientation, flipped when
/// the transcript is on '-' (the writer reverse-complements it). This, not
/// `align.is_reverse` (transcript strand vs read strand), is what the mate's
/// 0x20 flag must report.
fn output_reverse(eval: &ReadEval, m: &ExonChainMatch) -> bool {
    (eval.flags & FLAG_REVERSE != 0) ^ (m.align.strand == '-')
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
            let Some(read_match) = read.matches.get(&tid) else {
                continue;
            };
            let Some(mate_match) = mate.matches.get(&tid) else {
                continue;
            };

            let mate_info_for_read = MateInfo {
                tid,
                pos: align_pos(mate_match),
                is_reverse: output_reverse(mate, mate_match),
            };
            let mate_info_for_mate = MateInfo {
                tid,
                pos: align_pos(read_match),
                is_reverse: output_reverse(read, read_match),
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

        let Some(read_match) = read.matches.get(&r_tid) else {
            return groups;
        };
        let Some(mate_match) = mate.matches.get(&m_tid) else {
            return groups;
        };

        let mate_info_for_read = MateInfo {
            tid: m_tid,
            pos: align_pos(mate_match),
            is_reverse: output_reverse(mate, mate_match),
        };
        let mate_info_for_mate = MateInfo {
            tid: r_tid,
            pos: align_pos(read_match),
            is_reverse: output_reverse(read, read_match),
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
                    segs.push(alignment::Segment {
                        start: exon_start,
                        end: ref_pos,
                    });
                }
                ref_pos = ref_pos.saturating_add(len);
                exon_start = ref_pos;
            }
            _ => {}
        }
    }
    if ref_pos > exon_start {
        segs.push(alignment::Segment {
            start: exon_start,
            end: ref_pos,
        });
    }
    segs
}
