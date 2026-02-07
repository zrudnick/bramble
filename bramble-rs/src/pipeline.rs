use crate::alignment;
use crate::bam_input::BamInput;
use crate::cli::Args;
use crate::evaluate::{CigarOp, ExonChainMatch, ReadAln, ReadEvaluator};
use crate::g2t::G2TTree;
use crate::types::{ReadId, RefId, Tid};
use anyhow::Result;
use crossfire::mpmc;
use noodles::{bam, sam};
use noodles::core::Position;
use sam::alignment::io::Write as _;
use sam::alignment::record::cigar::{op::Kind as CigarKind, Op as SamCigarOp};
use sam::alignment::record::data::field::Tag;
use sam::alignment::record::Flags;
use sam::alignment::record::MappingQuality;
use sam::alignment::record_buf::{data::field::Value, Cigar as SamCigar, Data as SamData, QualityScores, Sequence};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::sync::Mutex;
use std::thread;

const UNORDERED_FLUSH_GROUPS: usize = 8;

const FR_STRAND: bool = true;
const RF_STRAND: bool = false;

#[derive(Debug, Default)]
pub struct Stats {
    pub total_reads: u64,
    pub unmapped_reads: u64,
    pub read_groups: u64,
    pub total_exons: u64,
}

#[derive(Debug)]
struct ReadRec {
    id: ReadId,
    record: bam::Record,
    segs: Vec<alignment::Segment>,
}

#[derive(Debug)]
struct ReadEval {
    record_idx: usize,
    matches: HashMap<Tid, ExonChainMatch>,
    read_len: usize,
    flags: Flags,
    alignment_start: Option<u32>,
    mate_alignment_start: Option<u32>,
    ref_id: Option<RefId>,
    mate_ref_id: Option<RefId>,
    hit_index: i32,
}

#[derive(Debug, Clone)]
struct MateInfo {
    tid: Tid,
    pos: u32,
    read_len: usize,
    is_reverse: bool,
}

#[derive(Debug, Clone)]
struct OutputEntry {
    record_idx: usize,
    tid: Tid,
    align: ExonChainMatch,
    read_len: usize,
    mate: Option<MateInfo>,
    is_first: bool,
    is_last: bool,
    same_transcript: bool,
    nh: u32,
}

pub fn run(
    args: &Args,
    out_header: &sam::Header,
    bam: &mut BamInput,
    g2t: &G2TTree,
    fasta: Option<&crate::fasta::FastaDb>,
) -> Result<Stats> {
    let out_file = File::create(&args.out_bam)?;
    let mut writer = bam::io::Writer::new(out_file);
    writer.write_header(out_header)?;

    let evaluator = ReadEvaluator {
        long_reads: args.long,
        use_fasta: fasta.is_some(),
    };

    let mut stats = Stats::default();
    let mut next_read_id: ReadId = 0;

    let mut current_name: Option<String> = None;
    let mut group: Vec<ReadRec> = Vec::new();

    if args.threads > 1 && args.unordered {
        crossfire::detect_backoff_cfg();
        let worker_count = args.threads as usize;
        let cap = worker_count.saturating_mul(4).max(8);
        let (tx_work, rx_work) = mpmc::bounded_blocking::<WorkItem>(cap);
        let writer_mutex = Mutex::new(writer);
        let mut group_idx: usize = 0;
        let mut total_groups: usize = 0;

        let evaluator_ref = &evaluator;
        let g2t_ref = g2t;
        let flush_records = args.unordered_flush_records.max(1);
        thread::scope(|scope| -> Result<()> {
            for _ in 0..worker_count {
                let rx_work = rx_work.clone();
                let evaluator_ref = evaluator_ref;
                let g2t_ref = g2t_ref;
                let writer_mutex = &writer_mutex;
                let flush_records = flush_records;
                let out_header = out_header;
                scope.spawn(move || {
                    let mut buffered: Vec<sam::alignment::RecordBuf> = Vec::new();
                    let mut buffered_groups: usize = 0;
                    while let Ok(item) = rx_work.recv() {
                        let result = process_group_records(&item.group, g2t_ref, evaluator_ref);
                        if let Ok(records) = result {
                            buffered_groups += 1;
                            buffered.extend(records);
                            if buffered_groups >= UNORDERED_FLUSH_GROUPS
                                || buffered.len() >= flush_records
                            {
                                if let Ok(mut writer) = writer_mutex.lock() {
                                    for record in buffered.drain(..) {
                                        let _ =
                                            writer.write_alignment_record(out_header, &record);
                                    }
                                } else {
                                    buffered.clear();
                                }
                                buffered_groups = 0;
                            }
                        }
                    }
                    if !buffered.is_empty() {
                        if let Ok(mut writer) = writer_mutex.lock() {
                            for record in buffered.drain(..) {
                                let _ = writer.write_alignment_record(out_header, &record);
                            }
                        }
                    }
                });
            }

            for result in bam.reader.records() {
                let record = result?;
                stats.total_reads += 1;
                if record.flags().is_unmapped() {
                    stats.unmapped_reads += 1;
                }

                let exons = alignment::extract_exons(&record)?;
                stats.total_exons += exons.len() as u64;

                let read_id = next_read_id;
                next_read_id = next_read_id.saturating_add(1);

                let name = record.name().map(|n| n.to_string()).unwrap_or_default();
                match &current_name {
                    None => {
                        current_name = Some(name);
                        group.push(ReadRec {
                            id: read_id,
                            record,
                            segs: exons,
                        });
                    }
                    Some(curr) if *curr == name => {
                        group.push(ReadRec {
                            id: read_id,
                            record,
                            segs: exons,
                        });
                    }
                    Some(_) => {
                        stats.read_groups += 1;
                        tx_work.send(WorkItem {
                            idx: group_idx,
                            group: std::mem::take(&mut group),
                        })?;
                        total_groups += 1;
                        group_idx += 1;
                        current_name = Some(name);
                        group.push(ReadRec {
                            id: read_id,
                            record,
                            segs: exons,
                        });
                    }
                }
            }

            if !group.is_empty() {
                stats.read_groups += 1;
                tx_work.send(WorkItem {
                    idx: group_idx,
                    group: std::mem::take(&mut group),
                })?;
                total_groups += 1;
                group_idx += 1;
            }

            drop(tx_work);

            Ok(())
        })?;

        return Ok(stats);
    }

    if args.threads > 1 {
        crossfire::detect_backoff_cfg();
        let worker_count = args.threads as usize;
        let cap = worker_count.saturating_mul(4).max(8);
        let (tx_work, rx_work) = mpmc::bounded_blocking::<WorkItem>(cap);
        let (tx_res, rx_res) = mpmc::unbounded_blocking::<ResultItem>();

        let mut group_idx: usize = 0;
        let mut total_groups: usize = 0;

        let evaluator_ref = &evaluator;
        let g2t_ref = g2t;
        thread::scope(|scope| -> Result<()> {
            for _ in 0..worker_count {
                let rx_work = rx_work.clone();
                let tx_res = tx_res.clone();
                let evaluator_ref = evaluator_ref;
                let g2t_ref = g2t_ref;
                scope.spawn(move || {
                    while let Ok(item) = rx_work.recv() {
                        let result = process_group_records(&item.group, g2t_ref, evaluator_ref);
                        let _ = tx_res.send(ResultItem { idx: item.idx, result });
                    }
                });
            }
            drop(tx_res);

            for result in bam.reader.records() {
                let record = result?;
                stats.total_reads += 1;
                if record.flags().is_unmapped() {
                    stats.unmapped_reads += 1;
                }

                let exons = alignment::extract_exons(&record)?;
                stats.total_exons += exons.len() as u64;

                let read_id = next_read_id;
                next_read_id = next_read_id.saturating_add(1);

                let name = record.name().map(|n| n.to_string()).unwrap_or_default();
                match &current_name {
                    None => {
                        current_name = Some(name);
                        group.push(ReadRec {
                            id: read_id,
                            record,
                            segs: exons,
                        });
                    }
                    Some(curr) if *curr == name => {
                        group.push(ReadRec {
                            id: read_id,
                            record,
                            segs: exons,
                        });
                    }
                    Some(_) => {
                        stats.read_groups += 1;
                        tx_work.send(WorkItem {
                            idx: group_idx,
                            group: std::mem::take(&mut group),
                        })?;
                        total_groups += 1;
                        group_idx += 1;
                        current_name = Some(name);
                        group.push(ReadRec {
                            id: read_id,
                            record,
                            segs: exons,
                        });
                    }
                }
            }

            if !group.is_empty() {
                stats.read_groups += 1;
                tx_work.send(WorkItem {
                    idx: group_idx,
                    group: std::mem::take(&mut group),
                })?;
                total_groups += 1;
                group_idx += 1;
            }

            drop(tx_work);

            let mut pending: BTreeMap<usize, Result<Vec<sam::alignment::RecordBuf>>> =
                BTreeMap::new();
            let mut next_idx = 0usize;
            let mut written = 0usize;

            while written < total_groups {
                let res = rx_res
                    .recv()
                    .map_err(|_| anyhow::anyhow!("worker result channel closed"))?;
                pending.insert(res.idx, res.result);
                while let Some(result) = pending.remove(&next_idx) {
                    let records = result?;
                    for record in records {
                        writer.write_alignment_record(out_header, &record)?;
                    }
                    next_idx += 1;
                    written += 1;
                }
            }

            Ok(())
        })?;

        return Ok(stats);
    }

    for result in bam.reader.records() {
        let record = result?;
        stats.total_reads += 1;
        if record.flags().is_unmapped() {
            stats.unmapped_reads += 1;
        }

        let exons = alignment::extract_exons(&record)?;
        stats.total_exons += exons.len() as u64;

        let read_id = next_read_id;
        next_read_id = next_read_id.saturating_add(1);

        let name = record.name().map(|n| n.to_string()).unwrap_or_default();
        match &current_name {
            None => {
                current_name = Some(name);
                group.push(ReadRec {
                    id: read_id,
                    record,
                    segs: exons,
                });
            }
            Some(curr) if *curr == name => {
                group.push(ReadRec {
                    id: read_id,
                    record,
                    segs: exons,
                });
            }
            Some(_) => {
                stats.read_groups += 1;
                let records = process_group_records(&group, g2t, &evaluator)?;
                for record in records {
                    writer.write_alignment_record(out_header, &record)?;
                }
                group.clear();
                current_name = Some(name);
                group.push(ReadRec {
                    id: read_id,
                    record,
                    segs: exons,
                });
            }
        }
    }

    if !group.is_empty() {
        stats.read_groups += 1;
        let records = process_group_records(&group, g2t, &evaluator)?;
        for record in records {
            writer.write_alignment_record(out_header, &record)?;
        }
    }

    Ok(stats)
}

struct WorkItem {
    idx: usize,
    group: Vec<ReadRec>,
}

struct ResultItem {
    idx: usize,
    result: Result<Vec<sam::alignment::RecordBuf>>,
}

fn process_group_records(
    group: &[ReadRec],
    g2t: &G2TTree,
    evaluator: &ReadEvaluator,
) -> Result<Vec<sam::alignment::RecordBuf>> {
    let mut output: Vec<sam::alignment::RecordBuf> = Vec::new();

    let mut read_evals: Vec<ReadEval> = Vec::new();
    let shared_seq: Option<Vec<u8>> = group
        .iter()
        .find(|rec| !rec.record.flags().is_unmapped())
        .map(|rec| rec.record.sequence().iter().collect());

    for (idx, rec) in group.iter().enumerate() {
        if rec.record.flags().is_unmapped() {
            continue;
        }

        let refid = match rec.record.reference_sequence_id() {
            Some(Ok(id)) => Some(id as RefId),
            Some(Err(_)) | None => None,
        };
        if refid.is_none() {
            continue;
        }

        let mate_refid = match rec.record.mate_reference_sequence_id() {
            Some(Ok(id)) => Some(id as RefId),
            Some(Err(_)) | None => None,
        };

        let (strand, strand_from_tag) = splice_strand(&rec.record);

        let read = ReadAln {
            strand,
            strand_from_tag,
            refid: refid.unwrap(),
            nh: 0,
            segs: rec.segs.clone(),
            record: rec.record.clone(),
        };

        let matches: HashMap<Tid, ExonChainMatch> = evaluator
            .evaluate(&read, rec.id, g2t, shared_seq.as_deref())
            .into_iter()
            .filter(|(_, m)| m.align.cigar.is_some())
            .collect();
        let read_len = rec.record.sequence().len();
        let hit_index = get_hit_index(&rec.record);

        read_evals.push(ReadEval {
            record_idx: idx,
            matches,
            read_len,
            flags: rec.record.flags(),
            alignment_start: get_alignment_start(&rec.record),
            mate_alignment_start: get_mate_start(&rec.record),
            ref_id: refid,
            mate_ref_id: mate_refid,
            hit_index,
        });
    }

    if read_evals.is_empty() {
        return Ok(output);
    }

    let pairs = find_mate_pairs(&read_evals);
    let mut paired = vec![false; read_evals.len()];

    let mut output_groups: Vec<Vec<OutputEntry>> = Vec::new();

    for (i, j) in pairs {
        paired[i] = true;
        paired[j] = true;

        let (r_idx, m_idx) = assign_pair_order(i, j, &read_evals);
        let read = &read_evals[r_idx];
        let mate = &read_evals[m_idx];

        output_groups.extend(build_paired_groups(read, mate));
    }

    for (idx, read) in read_evals.iter().enumerate() {
        if paired[idx] {
            continue;
        }
        output_groups.extend(build_unpaired_groups(read));
    }

    let new_nh = output_groups.len() as u32;
    let mut entries: Vec<OutputEntry> = output_groups.into_iter().flatten().collect();
    entries.sort_by(|a, b| {
        (
            a.record_idx,
            a.tid,
            a.is_first,
            a.is_last,
            a.align.align.hit_index,
        )
            .cmp(&(
                b.record_idx,
                b.tid,
                b.is_first,
                b.is_last,
                b.align.align.hit_index,
            ))
    });

    assign_hit_indices(&mut entries);
    for entry in entries {
        let is_read1 = entry.mate.is_none() || entry.is_first;
        let nh = if is_read1 { new_nh } else { entry.nh };
        let mapq = if is_read1 {
            get_mapq(nh, evaluator.long_reads)
        } else {
            0
        };
        let rec = &group[entry.record_idx].record;
        let out = build_projected_record(rec, &entry, nh, mapq, evaluator.long_reads)?;
        output.push(out);
    }

    Ok(output)
}

fn get_mapq(nh: u32, long_reads: bool) -> u8 {
    if !long_reads {
        if nh == 1 {
            255
        } else if nh == 2 {
            3
        } else if nh == 3 || nh == 4 {
            1
        } else {
            0
        }
    } else if nh > 1 {
        0
    } else {
        3
    }
}

fn assign_hit_indices(entries: &mut [OutputEntry]) {
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
            } else if (score - best_score).abs() <= 1e-12 {
                if let Some(current) = best_idx {
                    if entry.tid < entries[current].tid {
                        best_idx = Some(entry_idx);
                    }
                }
            }
        }

        let best_idx = best_idx.unwrap_or(indices[0]);
        for entry_idx in indices.iter().copied() {
            entries[entry_idx].align.align.primary_alignment = entry_idx == best_idx;
        }
    }
}

fn get_alignment_start(record: &bam::Record) -> Option<u32> {
    record
        .alignment_start()
        .and_then(|res| res.ok())
        .map(|pos| pos.get() as u32)
}

fn splice_strand(record: &bam::Record) -> (char, bool) {
    let xs_tag = Tag::new(b'X', b'S');
    let ts_tag = Tag::new(b't', b's');

    if let Some(c) = get_char_tag(record, xs_tag) {
        let strand = c as char;
        if strand == '+' || strand == '-' {
            return (strand, true);
        }
    }

    if let Some(c) = get_char_tag(record, ts_tag) {
        let mut strand = c as char;
        if strand == '+' || strand == '-' {
            if record.flags().is_reverse_complemented() {
                strand = if strand == '+' { '-' } else { '+' };
            }
            return (strand, true);
        }
    }
    // Match C++ behavior: if no splice-strand tags, infer from library layout.
    let flags = record.flags();
    let rev = flags.is_reverse_complemented();
    if FR_STRAND || RF_STRAND {
        let strand = if flags.is_segmented() {
            if flags.is_first_segment() {
                if (RF_STRAND && rev) || (FR_STRAND && !rev) {
                    '+'
                } else {
                    '-'
                }
            } else {
                if (RF_STRAND && rev) || (FR_STRAND && !rev) {
                    '-'
                } else {
                    '+'
                }
            }
        } else if (RF_STRAND && rev) || (FR_STRAND && !rev) {
            '+'
        } else {
            '-'
        };
        (strand, false)
    } else if rev {
        ('-', false)
    } else {
        ('+', false)
    }
}

fn get_char_tag(record: &bam::Record, tag: Tag) -> Option<u8> {
    let data = record.data();
    let value = data.get(&tag)?;
    let value = value.ok()?;
    match value {
        noodles::sam::alignment::record::data::field::Value::Character(c) => Some(c),
        noodles::sam::alignment::record::data::field::Value::String(s) => s.first().copied(),
        _ => None,
    }
}

fn get_mate_start(record: &bam::Record) -> Option<u32> {
    record
        .mate_alignment_start()
        .and_then(|res| res.ok())
        .map(|pos| pos.get() as u32)
}

fn get_hit_index(record: &bam::Record) -> i32 {
    let data = record.data();
    let value = data.get(&Tag::HIT_INDEX);
    if let Some(Ok(v)) = value {
        v.as_int().unwrap_or(0) as i32
    } else {
        0
    }
}

fn find_mate_pairs(reads: &[ReadEval]) -> Vec<(usize, usize)> {
    let mut pending: HashMap<(u32, i32), usize> = HashMap::new();
    let mut pairs = Vec::new();

    for (idx, read) in reads.iter().enumerate() {
        if !read.flags.is_segmented() || read.flags.is_mate_unmapped() {
            continue;
        }

        let (Some(read_start), Some(mate_start), Some(ref_id), Some(mate_ref_id)) = (
            read.alignment_start,
            read.mate_alignment_start,
            read.ref_id,
            read.mate_ref_id,
        ) else {
            continue;
        };

        if ref_id != mate_ref_id {
            continue;
        }

        let key = (read_start, read.hit_index);
        if mate_start <= read_start {
            let mate_key = (mate_start, read.hit_index);
            if let Some(mate_idx) = pending.remove(&mate_key) {
                pairs.push((idx, mate_idx));
            }
        } else {
            pending.entry(key).or_insert(idx);
        }
    }

    pairs
}

fn assign_pair_order(i: usize, j: usize, reads: &[ReadEval]) -> (usize, usize) {
    let flags_i = reads[i].flags;
    let flags_j = reads[j].flags;

    let i_first = flags_i.is_first_segment();
    let i_last = flags_i.is_last_segment();
    let j_first = flags_j.is_first_segment();
    let j_last = flags_j.is_last_segment();

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

fn build_paired_groups(read: &ReadEval, mate: &ReadEval) -> Vec<Vec<OutputEntry>> {
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
                read_len: mate.read_len,
                is_reverse: mate_match.align.is_reverse,
            };
            let mate_info_for_mate = MateInfo {
                tid,
                pos: align_pos(read_match),
                read_len: read.read_len,
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
            read_len: mate.read_len,
            is_reverse: mate_match.align.is_reverse,
        };
        let mate_info_for_mate = MateInfo {
            tid: r_tid,
            pos: align_pos(read_match),
            read_len: read.read_len,
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

fn build_unpaired_groups(read: &ReadEval) -> Vec<Vec<OutputEntry>> {
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

fn align_pos(match_info: &ExonChainMatch) -> u32 {
    if match_info.align.strand == '+' {
        match_info.align.fwpos
    } else {
        match_info.align.rcpos
    }
}

fn compute_template_length(
    pos: u32,
    mate_pos: u32,
    read_len: usize,
    mate_len: usize,
    is_first: bool,
    same_transcript: bool,
) -> i32 {
    if !same_transcript {
        return 0;
    }

    let pos = pos as i32;
    let mate_pos = mate_pos as i32;
    let read_len = read_len as i32;
    let mate_len = mate_len as i32;

    if is_first {
        if mate_pos > pos {
            (mate_pos + mate_len) - pos
        } else {
            -((pos + read_len) - mate_pos)
        }
    } else if pos > mate_pos {
        -((pos + mate_len) - mate_pos)
    } else {
        (mate_pos + read_len) - pos
    }
}

fn build_projected_record(
    record: &bam::Record,
    entry: &OutputEntry,
    nh: u32,
    mapq: u8,
    long_reads: bool,
) -> Result<sam::alignment::RecordBuf> {
    let mut out = sam::alignment::RecordBuf::default();

    if let Some(name) = record.name() {
        *out.name_mut() = Some(name.to_vec().into());
    }

    let original_reverse = record.flags().is_reverse_complemented();
    let align_reverse = entry.align.align.is_reverse;
    let has_sequence = !record.sequence().is_empty();
    let reverse_action = if long_reads {
        !original_reverse && align_reverse
    } else {
        original_reverse && has_sequence
    };
    let output_reverse = if reverse_action {
        !original_reverse
    } else {
        original_reverse
    };

    let mut flags = record.flags();
    flags.remove(
        Flags::UNMAPPED
            | Flags::MATE_UNMAPPED
            | Flags::MATE_REVERSE_COMPLEMENTED
            | Flags::PROPERLY_SEGMENTED
            | Flags::SEGMENTED
            | Flags::FIRST_SEGMENT
            | Flags::LAST_SEGMENT,
    );

    if entry.align.align.primary_alignment {
        flags.remove(Flags::SECONDARY);
    } else {
        flags.insert(Flags::SECONDARY);
    }

    if output_reverse {
        flags.insert(Flags::REVERSE_COMPLEMENTED);
    } else {
        flags.remove(Flags::REVERSE_COMPLEMENTED);
    }

    let pos0 = if entry.align.align.strand == '+' {
        entry.align.align.fwpos
    } else {
        entry.align.align.rcpos
    };

    let pos1 = pos0.saturating_add(1) as usize;
    let alignment_start = Position::try_from(pos1)
        .map_err(|_| anyhow::anyhow!("alignment start out of range: {pos1}"))?;

    *out.reference_sequence_id_mut() = Some(entry.tid as usize);
    *out.alignment_start_mut() = Some(alignment_start);

    let mapping_quality = if mapq == 255 {
        None
    } else {
        Some(MappingQuality::try_from(mapq).unwrap_or(MappingQuality::MIN))
    };
    *out.mapping_quality_mut() = mapping_quality;

    let (mut cigar, nm) = update_cigar(record, entry.align.align.cigar.as_ref().unwrap())?;
    if reverse_action {
        cigar.as_mut().reverse();
    }
    *out.cigar_mut() = cigar;

    if let Some(mate) = entry.mate.as_ref() {
        flags.insert(Flags::SEGMENTED);
        flags.insert(Flags::PROPERLY_SEGMENTED);
        if entry.is_first {
            flags.insert(Flags::FIRST_SEGMENT);
        }
        if entry.is_last {
            flags.insert(Flags::LAST_SEGMENT);
        }
        if mate.is_reverse {
            flags.insert(Flags::MATE_REVERSE_COMPLEMENTED);
        } else {
            flags.remove(Flags::MATE_REVERSE_COMPLEMENTED);
        }

        let mate_pos = mate.pos.saturating_add(1) as usize;
        let mate_alignment_start = Position::try_from(mate_pos)
            .map_err(|_| anyhow::anyhow!("mate alignment start out of range: {mate_pos}"))?;

        *out.mate_reference_sequence_id_mut() = Some(mate.tid as usize);
        *out.mate_alignment_start_mut() = Some(mate_alignment_start);
        *out.template_length_mut() = compute_template_length(
            pos0,
            mate.pos,
            entry.read_len,
            mate.read_len,
            entry.is_first,
            entry.same_transcript,
        );
    } else {
        flags.remove(
            Flags::SEGMENTED
                | Flags::PROPERLY_SEGMENTED
                | Flags::MATE_UNMAPPED
                | Flags::MATE_REVERSE_COMPLEMENTED
                | Flags::FIRST_SEGMENT
                | Flags::LAST_SEGMENT,
        );
        *out.mate_reference_sequence_id_mut() = None;
        *out.mate_alignment_start_mut() = None;
        *out.template_length_mut() = 0;
    }

    *out.flags_mut() = flags;

    let mut seq: Vec<u8> = record.sequence().iter().collect();
    let mut qual: Vec<u8> = record.quality_scores().iter().collect();

    if reverse_action {
        reverse_complement(&mut seq);
        qual.reverse();
    }

    *out.sequence_mut() = Sequence::from(seq);
    *out.quality_scores_mut() = QualityScores::from(qual);

    let mut data: SamData = record
        .data()
        .try_into()
        .map_err(|e: std::io::Error| anyhow::anyhow!(e))?;
    remove_unwanted_tags(&mut data);
    data.insert(Tag::ALIGNMENT_HIT_COUNT, Value::from(nh as i32));
    data.insert(Tag::HIT_INDEX, Value::from(entry.align.align.hit_index));
    data.insert(
        Tag::new(b'X', b'S'),
        Value::Character(entry.align.align.strand as u8),
    );

    let gn_as = extract_alignment_score(record).unwrap_or(0) as f64;
    let score = if !long_reads {
        (1.0 + entry.align.align.similarity_score).powi(3) * 100.0
    } else {
        (gn_as + (entry.align.align.clip_score as f64)) * entry.align.align.similarity_score
    };
    let as_tag = score as i32;
    data.insert(Tag::ALIGNMENT_SCORE, Value::from(as_tag));
    data.insert(Tag::EDIT_DISTANCE, Value::from(nm));

    *out.data_mut() = data;

    Ok(out)
}


fn extract_alignment_score(record: &bam::Record) -> Option<i32> {
    let data = record.data();
    let value = data.get(&Tag::ALIGNMENT_SCORE)?;
    let value = value.ok()?;
    value.as_int().map(|v| v as i32)
}

fn reverse_complement(seq: &mut [u8]) {
    seq.reverse();
    for base in seq.iter_mut() {
        *base = match *base {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'N' | b'n' => b'N',
            _ => b'N',
        };
    }
}

fn remove_unwanted_tags(data: &mut SamData) {
    let tags = [
        Tag::ALIGNMENT_HIT_COUNT,
        Tag::HIT_INDEX,
        Tag::ALIGNMENT_SCORE,
        Tag::EDIT_DISTANCE,
        Tag::CIGAR,
        Tag::new(b'X', b'S'),
        Tag::new(b'S', b'A'),
        Tag::new(b'm', b's'),
        Tag::new(b'n', b'n'),
        Tag::new(b't', b's'),
        Tag::new(b't', b'p'),
        Tag::new(b'c', b'm'),
        Tag::new(b's', b'1'),
        Tag::new(b's', b'2'),
        Tag::new(b'd', b'e'),
        Tag::new(b'r', b'l'),
    ];

    for tag in tags {
        data.remove(&tag);
    }
}

fn update_cigar(
    record: &bam::Record,
    ideal: &crate::evaluate::Cigar,
) -> Result<(SamCigar, i32)> {
    let mut real_ops: Vec<SamCigarOp> = Vec::new();
    for result in record.cigar().iter() {
        real_ops.push(result?);
    }

    Ok(update_cigar_ops(&real_ops, ideal))
}

fn update_cigar_ops(real_ops: &[SamCigarOp], ideal: &crate::evaluate::Cigar) -> (SamCigar, i32) {
    let (front_hard, front_soft) = front_clip_lengths(real_ops);

    let mut real_expanded = expand_real_cigar(real_ops);
    let ideal_runs = ideal_runs(ideal);
    let merged_ideal = merge_indels(&ideal_runs);
    let mut ideal_expanded = expand_runs(&merged_ideal);

    pad_ideal_for_leading_clips(&mut ideal_expanded, front_hard, front_soft);
    align_expanded_cigars(&mut real_expanded, &mut ideal_expanded);

    let max_len = real_expanded.len().max(ideal_expanded.len());
    let mut merged: Vec<char> = Vec::with_capacity(max_len);
    for i in 0..max_len {
        let real_op = real_expanded.get(i).copied().unwrap_or('*');
        let ideal_op = ideal_expanded.get(i).copied().unwrap_or('*');
        merged.push(merge_ops(real_op, ideal_op));
    }

    let mut nm = 0i32;
    let compressed = compress_cigar(&merged, &mut nm);
    let final_runs = merge_indels(&compressed);
    let cigar = runs_to_sam_cigar(&final_runs);

    (cigar, nm)
}

#[allow(dead_code)]
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

fn expand_real_cigar(ops: &[SamCigarOp]) -> Vec<char> {
    let mut expanded = Vec::new();
    for op in ops {
        if op.kind() == CigarKind::Skip {
            continue;
        }
        let ch = cigar_kind_to_char(op.kind());
        expanded.extend(std::iter::repeat_n(ch, op.len()));
    }
    expanded
}

fn cigar_kind_to_char(kind: CigarKind) -> char {
    match kind {
        CigarKind::Match => 'M',
        CigarKind::Insertion => 'I',
        CigarKind::Deletion => 'D',
        CigarKind::Skip => 'N',
        CigarKind::SoftClip => 'S',
        CigarKind::HardClip => 'H',
        CigarKind::Pad => 'P',
        CigarKind::SequenceMatch => '=',
        CigarKind::SequenceMismatch => 'X',
    }
}

fn ideal_runs(cigar: &crate::evaluate::Cigar) -> Vec<(u32, char)> {
    let mut runs = Vec::new();
    for (len, op) in cigar.ops.iter() {
        if *len == 0 {
            continue;
        }
        let ch = match op {
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
        };
        runs.push((*len, ch));
    }
    runs
}

fn merge_indels(runs: &[(u32, char)]) -> Vec<(u32, char)> {
    let mut result: Vec<(u32, char)> = Vec::new();
    let mut i_count: u32 = 0;
    let mut d_count: u32 = 0;

    let push_run = |out: &mut Vec<(u32, char)>, len: u32, op: char| {
        if len == 0 {
            return;
        }
        if let Some(last) = out.last_mut()
            && last.1 == op
        {
            last.0 = last.0.saturating_add(len);
        } else {
            out.push((len, op));
        }
    };

    let flush = |out: &mut Vec<(u32, char)>, i_count: &mut u32, d_count: &mut u32| {
        if *i_count == 0 && *d_count == 0 {
            return;
        }
        let overlap = (*i_count).min(*d_count);
        if overlap > 0 {
            push_run(out, overlap, 'M');
            *i_count -= overlap;
            *d_count -= overlap;
        }
        if *i_count > 0 {
            push_run(out, *i_count, 'I');
        }
        if *d_count > 0 {
            push_run(out, *d_count, 'D');
        }
        *i_count = 0;
        *d_count = 0;
    };

    for (len, op) in runs.iter().copied() {
        if op == 'I' {
            i_count += len;
        } else if op == 'D' {
            d_count += len;
        } else {
            flush(&mut result, &mut i_count, &mut d_count);
            push_run(&mut result, len, op);
        }
    }
    flush(&mut result, &mut i_count, &mut d_count);

    result
}

fn expand_runs(runs: &[(u32, char)]) -> Vec<char> {
    let mut expanded = Vec::new();
    for (len, op) in runs {
        expanded.extend(std::iter::repeat_n(*op, *len as usize));
    }
    expanded
}

fn pad_ideal_for_leading_clips(ideal: &mut Vec<char>, front_hard: u32, front_soft: u32) {
    let total_padding = front_hard.saturating_add(front_soft);
    if total_padding == 0 {
        return;
    }

    let mut padded = Vec::with_capacity(total_padding as usize + ideal.len());

    padded.extend(std::iter::repeat_n('_', front_hard as usize));

    for i in front_hard..front_soft {
        let idx = i as usize;
        if idx < ideal.len() && matches!(ideal[idx], ',' | '.' | '/' | ';') {
            continue;
        }
        padded.push('_');
    }

    padded.extend_from_slice(ideal);
    *ideal = padded;
}

fn align_expanded_cigars(real: &mut Vec<char>, ideal: &mut Vec<char>) {
    let mut aligned_real = Vec::new();
    let mut aligned_ideal = Vec::new();

    let mut real_pos = 0usize;
    let mut ideal_pos = 0usize;
    let max_iterations = (real.len() + ideal.len()) * 2;
    let mut iterations = 0usize;

    while real_pos < real.len() || ideal_pos < ideal.len() {
        iterations += 1;
        if iterations > max_iterations {
            break;
        }

        if real_pos >= real.len() {
            aligned_real.push('_');
            aligned_ideal.push(ideal[ideal_pos]);
            ideal_pos += 1;
            continue;
        }
        if ideal_pos >= ideal.len() {
            aligned_real.push(real[real_pos]);
            aligned_ideal.push('_');
            real_pos += 1;
            continue;
        }

        let r = real[real_pos];
        let i = ideal[ideal_pos];

        if i == '.' {
            aligned_real.push('_');
            aligned_ideal.push(i);
            ideal_pos += 1;
        } else if r == 'I' {
            aligned_real.push(r);
            aligned_ideal.push('_');
            real_pos += 1;
        } else if i == 'D' {
            aligned_real.push('_');
            aligned_ideal.push(i);
            ideal_pos += 1;
        } else {
            aligned_real.push(r);
            aligned_ideal.push(i);
            real_pos += 1;
            ideal_pos += 1;
        }
    }

    *real = aligned_real;
    *ideal = aligned_ideal;
}

fn merge_ops(real_op: char, ideal_op: char) -> char {
    if real_op == 'I' && ideal_op == '_' {
        return 'I';
    }

    if (real_op == 'M' || real_op == 'S') && ideal_op == ';' {
        return 'S';
    }
    if (real_op == 'M' || real_op == 'S') && ideal_op == ',' {
        return 'M';
    }
    if (real_op == 'M' || real_op == 'S') && ideal_op == '/' {
        return 'I';
    }
    if (real_op == 'M' || real_op == 'S') && ideal_op == '.' {
        return 'D';
    }

    if real_op == 'D' && ideal_op == ';' {
        return '_';
    }
    if real_op == 'D' && ideal_op == ',' {
        return 'D';
    }
    if real_op == 'D' && ideal_op == '/' {
        return '_';
    }
    if real_op == 'D' && ideal_op == '.' {
        return '_';
    }

    if real_op == 'I' && ideal_op == ';' {
        return 'S';
    }
    if real_op == 'I' && ideal_op == ',' {
        return 'I';
    }
    if real_op == 'D' && ideal_op == '/' {
        return '_';
    }
    if real_op == 'I' && ideal_op == '.' {
        return '_';
    }

    if ideal_op == ';' {
        return 'S';
    }
    if ideal_op == ',' {
        return 'M';
    }
    if ideal_op == '/' {
        return 'I';
    }
    if ideal_op == '.' {
        return 'D';
    }

    if ideal_op == '*' {
        return real_op;
    }
    if real_op == '*' {
        return ideal_op;
    }

    if real_op == 'H' {
        return 'H';
    }

    if real_op == 'D' && ideal_op == 'S' {
        return '_';
    }
    if real_op == 'I' && ideal_op == 'S' {
        return 'S';
    }
    if real_op == 'D' && ideal_op == 'I' {
        return '_';
    }

    if ideal_op == 'S' || ideal_op == 'D' || ideal_op == 'I' {
        return ideal_op;
    }

    if real_op == 'S' || real_op == 'D' || real_op == 'I' {
        return real_op;
    }

    if ideal_op == 'M' || ideal_op == '=' || ideal_op == 'X' {
        return 'M';
    }
    if real_op == 'M' || real_op == '=' || real_op == 'X' {
        return 'M';
    }

    if real_op == '_' {
        return ideal_op;
    }
    if ideal_op == '_' {
        return real_op;
    }

    if real_op != '*' {
        real_op
    } else {
        ideal_op
    }
}

fn compress_cigar(expanded: &[char], nm: &mut i32) -> Vec<(u32, char)> {
    let mut cleaned: Vec<char> = expanded.to_vec();
    if cleaned.len() >= 3 {
        for i in 1..(cleaned.len() - 1) {
            if cleaned[i] == 'I' {
                let mut has_s_before = false;
                let mut j = i;
                while j > 0 {
                    let c = cleaned[j - 1];
                    if c == 'S' {
                        has_s_before = true;
                        break;
                    }
                    if c != 'I' {
                        break;
                    }
                    j -= 1;
                }

                let mut has_s_after = false;
                let mut j = i;
                while j + 1 < cleaned.len() {
                    let c = cleaned[j + 1];
                    if c == 'S' {
                        has_s_after = true;
                        break;
                    }
                    if c != 'I' {
                        break;
                    }
                    j += 1;
                }

                if has_s_before && has_s_after {
                    cleaned[i] = 'S';
                }
            }
        }
    }

    let mut runs = Vec::new();
    let mut i = 0usize;
    while i < cleaned.len() {
        if cleaned[i] == '_' {
            i += 1;
            continue;
        }

        let op_char = cleaned[i];
        let mut run_len = 0u32;
        while i < cleaned.len() && (cleaned[i] == op_char || cleaned[i] == '_') {
            if cleaned[i] == op_char {
                run_len += 1;
            }
            i += 1;
        }

        if run_len == 0 {
            continue;
        }

        let out_char = match op_char {
            'M' | '=' | 'X' => 'M',
            'I' => {
                *nm += run_len as i32;
                'I'
            }
            'D' => {
                *nm += run_len as i32;
                'D'
            }
            'S' => 'S',
            'H' => 'H',
            'N' => 'N',
            'P' => 'P',
            _ => continue,
        };

        runs.push((run_len, out_char));
    }

    runs
}

fn runs_to_sam_cigar(runs: &[(u32, char)]) -> SamCigar {
    let ops: Vec<SamCigarOp> = runs
        .iter()
        .filter_map(|(len, op)| {
            if *len == 0 {
                return None;
            }
            let kind = match op {
                'M' => CigarKind::Match,
                'I' => CigarKind::Insertion,
                'D' => CigarKind::Deletion,
                'N' => CigarKind::Skip,
                'S' => CigarKind::SoftClip,
                'H' => CigarKind::HardClip,
                'P' => CigarKind::Pad,
                '=' => CigarKind::SequenceMatch,
                'X' => CigarKind::SequenceMismatch,
                _ => return None,
            };
            Some(SamCigarOp::new(kind, *len as usize))
        })
        .collect();

    ops.into_iter().collect()
}
