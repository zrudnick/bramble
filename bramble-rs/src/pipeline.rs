// pipeline.rs is shared between the library and binary targets.  Functions used only
// by one target look "unused" to the other — suppress those false-positive warnings.
#![allow(dead_code)]
use crate::alignment;
use crate::bam_input::BamInput;
use crate::cli::Args;
use crate::evaluate::{CigarOp, EvalContext, ExonChainMatch, ReadAln, ReadEvaluator};
use crate::g2t::G2TTree;
use crate::types::{HashMap, HashMapExt, HashSet, ReadId, RefId, Tid};
use anyhow::Result;
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use noodles::sam::alignment::record::cigar::{op::Kind as CigarKind, Op as SamCigarOp};
use noodles::sam::alignment::record_buf::Cigar as SamCigar;
use rust_htslib::bam as hts_bam;
use rust_htslib::bam::Read as HtsRead;
use rust_htslib::bam::record::{Aux, CigarString};
use std::collections::BTreeMap;
use std::sync::Mutex;
use std::thread;

const UNORDERED_FLUSH_GROUPS: usize = 8;
const PROGRESS_UPDATE_INTERVAL: u64 = 1000;
const BATCH_SIZE: usize = 64;

/// Tags to strip from output records (they are invalid in projected coordinates
/// or will be replaced with new values).
const SKIP_TAGS: &[[u8; 2]] = &[
    *b"NH", *b"HI", *b"AS", *b"NM", *b"CG",
    *b"XS", *b"SA", *b"ms", *b"nn", *b"ts",
    *b"tp", *b"cm", *b"s1", *b"s2", *b"de", *b"rl",
];

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
    record: hts_bam::Record,
    segs: Vec<alignment::Segment>,
    decoded_seq: Vec<u8>,
}

#[derive(Debug)]
pub(crate) struct ReadEval {
    pub(crate) record_idx: usize,
    pub(crate) matches: HashMap<Tid, ExonChainMatch>,
    pub(crate) read_len: usize,
    pub(crate) flags: u16,
    pub(crate) alignment_start: Option<u32>,
    pub(crate) mate_alignment_start: Option<u32>,
    pub(crate) ref_id: Option<RefId>,
    pub(crate) mate_ref_id: Option<RefId>,
    pub(crate) hit_index: i32,
}

// BAM flag constants
const FLAG_PAIRED: u16       = 0x1;
const FLAG_PROPER_PAIR: u16  = 0x2;
const FLAG_UNMAPPED: u16     = 0x4;
const FLAG_MATE_UNMAPPED: u16 = 0x8;
const FLAG_REVERSE: u16      = 0x10;
const FLAG_MATE_REVERSE: u16 = 0x20;
const FLAG_READ1: u16        = 0x40;
const FLAG_READ2: u16        = 0x80;
const FLAG_SECONDARY: u16    = 0x100;

#[derive(Debug, Clone)]
pub(crate) struct MateInfo {
    pub(crate) tid: Tid,
    pub(crate) pos: u32,
    pub(crate) read_len: usize,
    pub(crate) is_reverse: bool,
}

#[derive(Debug, Clone)]
pub(crate) struct OutputEntry {
    pub(crate) record_idx: usize,
    pub(crate) tid: Tid,
    pub(crate) align: ExonChainMatch,
    pub(crate) read_len: usize,
    pub(crate) mate: Option<MateInfo>,
    pub(crate) is_first: bool,
    pub(crate) is_last: bool,
    pub(crate) same_transcript: bool,
    pub(crate) nh: u32,
}

pub fn run(
    args: &Args,
    hts_header: &rust_htslib::bam::Header,
    bam: &mut BamInput,
    g2t: &G2TTree,
    fasta: Option<&crate::fasta::FastaDb>,
) -> Result<Stats> {
    let mut writer = hts_bam::Writer::from_path(&args.out_bam, hts_header, hts_bam::Format::Bam)?;
    writer.set_threads(4)?;

    let progress = if !args.quiet {
        let pb = ProgressBar::new_spinner();
        pb.set_draw_target(ProgressDrawTarget::stderr_with_hz(2));
        pb.set_style(
            ProgressStyle::default_spinner()
                .template("{spinner:.green} [{elapsed_precise}] {msg}")
                .expect("Failed to set progress bar template"),
        );
        pb.set_message("Processing alignments...");
        Some(pb)
    } else {
        None
    };

    let evaluator = ReadEvaluator { long_reads: args.long, use_fasta: fasta.is_some() };
    let fr = args.fr;
    let rf = args.rf;
    let paired_end = args.paired_end;
    let worker_count = args.threads as usize;
    let cap = worker_count.saturating_mul(4).max(8);
    let flush_records = args.unordered_flush_records.max(1);
    let evaluator_ref = &evaluator;
    let g2t_ref = g2t;

    // Always spawn a dedicated reader thread so BAM parsing overlaps with projection.
    // The reader thread reads, groups, and sends WorkItems; the main thread (and any
    // additional worker threads) handles projection and writing.

    let stats = if args.unordered {
        // Unordered: N worker threads write via a mutex; output order is not guaranteed.
        let (tx_work, rx_work) = flume::bounded::<WorkItem>(cap);
        let writer_mutex = Mutex::new(writer);

        thread::scope(|scope| -> Result<Stats> {
            for _ in 0..worker_count {
                let rx_work = rx_work.clone();
                let writer_mutex = &writer_mutex;
                scope.spawn(move || {
                    let mut ctx = EvalContext::new();
                    let mut buffered: Vec<hts_bam::Record> = Vec::new();
                    let mut buffered_groups: usize = 0;
                    while let Ok(item) = rx_work.recv() {
                        for (_, group) in &item.groups {
                            let result = process_group_records(
                                group, g2t_ref, evaluator_ref, fr, rf, paired_end, &mut ctx,
                            );
                            if let Ok(records) = result {
                                buffered_groups += 1;
                                buffered.extend(records);
                            }
                        }
                        if buffered_groups >= UNORDERED_FLUSH_GROUPS
                            || buffered.len() >= flush_records
                        {
                            if let Ok(mut w) = writer_mutex.lock() {
                                for record in buffered.drain(..) {
                                    let _ = w.write(&record);
                                }
                            } else {
                                buffered.clear();
                            }
                            buffered_groups = 0;
                        }
                    }
                    if !buffered.is_empty()
                        && let Ok(mut w) = writer_mutex.lock()
                    {
                        for record in buffered.drain(..) {
                            let _ = w.write(&record);
                        }
                    }
                });
            }

            // Reader thread: reads BAM, groups by name, sends WorkItems.
            let reader_jh =
                scope.spawn(|| read_and_group(&mut bam.reader, tx_work, &progress));

            // Main thread waits; workers and reader joined by scope exit.
            let stats = reader_jh
                .join()
                .map_err(|_| anyhow::anyhow!("reader thread panicked"))?;
            Ok(stats)
        })?
    } else if worker_count == 1 {
        // Single processing thread (main) + dedicated reader thread.
        // Main processes groups serially and writes in order as they arrive.
        let (tx_work, rx_work) = flume::bounded::<WorkItem>(cap);

        thread::scope(|scope| -> Result<Stats> {
            let reader_jh =
                scope.spawn(|| read_and_group(&mut bam.reader, tx_work, &progress));

            let mut ctx = EvalContext::new();
            while let Ok(item) = rx_work.recv() {
                for (_, group) in &item.groups {
                    let records = process_group_records(
                        group, g2t_ref, evaluator_ref, fr, rf, paired_end, &mut ctx,
                    )?;
                    for record in &records {
                        writer.write(record)?;
                    }
                }
            }

            let stats = reader_jh
                .join()
                .map_err(|_| anyhow::anyhow!("reader thread panicked"))?;
            Ok(stats)
        })?
    } else {
        // N worker threads + dedicated reader thread + main writes in order.
        let (tx_work, rx_work) = flume::bounded::<WorkItem>(cap);
        let (tx_res, rx_res) = flume::unbounded::<ResultItem>();

        thread::scope(|scope| -> Result<Stats> {
            for _ in 0..worker_count {
                let rx_work = rx_work.clone();
                let tx_res = tx_res.clone();
                scope.spawn(move || {
                    let mut ctx = EvalContext::new();
                    while let Ok(item) = rx_work.recv() {
                        let mut results: Vec<(usize, Result<Vec<hts_bam::Record>>)> = Vec::with_capacity(item.groups.len());
                        for (idx, group) in &item.groups {
                            let result = process_group_records(
                                group, g2t_ref, evaluator_ref, fr, rf, paired_end, &mut ctx,
                            );
                            results.push((*idx, result));
                        }
                        let _ = tx_res.send(ResultItem { results });
                    }
                });
            }
            drop(tx_res); // workers own the only remaining senders

            // Reader thread: send WorkItems; dropping tx_work when done signals workers.
            let reader_jh =
                scope.spawn(|| read_and_group(&mut bam.reader, tx_work, &progress));

            // Main: receive results and write in index order.
            let mut pending: BTreeMap<usize, Result<Vec<hts_bam::Record>>> =
                BTreeMap::new();
            let mut next_idx = 0usize;
            while let Ok(res) = rx_res.recv() {
                for (idx, result) in res.results {
                    pending.insert(idx, result);
                }
                while let Some(result) = pending.remove(&next_idx) {
                    let records = result?;
                    for record in &records {
                        writer.write(record)?;
                    }
                    next_idx += 1;
                }
            }

            let stats = reader_jh
                .join()
                .map_err(|_| anyhow::anyhow!("reader thread panicked"))?;
            Ok(stats)
        })?
    };

    if let Some(pb) = progress {
        pb.finish_with_message(format!("Completed: {} reads processed", stats.total_reads));
    }

    Ok(stats)
}

/// Reads BAM records, groups them by query name, batches groups, and sends
/// batches to `tx` as [`WorkItem`]s.  Returns accumulated statistics.
///
/// Dropping `tx` (when this function returns) closes the channel, signalling
/// downstream consumers to exit their recv loops.
fn read_and_group(
    reader: &mut hts_bam::Reader,
    tx: flume::Sender<WorkItem>,
    progress: &Option<ProgressBar>,
) -> Stats {
    let mut stats = Stats::default();
    let mut current_qname: Option<Vec<u8>> = None;
    let mut group: Vec<ReadRec> = Vec::new();
    let mut group_idx: usize = 0;
    let mut next_read_id: ReadId = 0;
    let mut batch: Vec<(usize, Vec<ReadRec>)> = Vec::with_capacity(BATCH_SIZE);

    let mut record = hts_bam::Record::new();
    loop {
        match reader.read(&mut record) {
            None => break,
            Some(Err(_)) => break,
            Some(Ok(())) => {}
        }
        stats.total_reads += 1;
        if record.is_unmapped() {
            stats.unmapped_reads += 1;
        }
        if let Some(pb) = progress
            && stats.total_reads.is_multiple_of(PROGRESS_UPDATE_INTERVAL)
        {
            pb.set_message(format!("Processed {} reads", stats.total_reads));
            pb.tick();
        }
        let exons = match alignment::extract_exons(&record) {
            Ok(e) => e,
            Err(_) => continue,
        };
        stats.total_exons += exons.len() as u64;
        let read_id = next_read_id;
        next_read_id = next_read_id.saturating_add(1);
        let qname = record.qname();
        let decoded_seq = record.seq().as_bytes();
        let same_name = current_qname.as_deref() == Some(qname);
        if same_name || current_qname.is_none() {
            if current_qname.is_none() {
                current_qname = Some(qname.to_vec());
            }
            group.push(ReadRec { id: read_id, record: record.clone(), segs: exons, decoded_seq });
        } else {
            // Name changed: flush current group into batch
            stats.read_groups += 1;
            batch.push((group_idx, std::mem::take(&mut group)));
            group_idx += 1;
            if batch.len() >= BATCH_SIZE
                && tx.send(WorkItem { groups: std::mem::replace(&mut batch, Vec::with_capacity(BATCH_SIZE)) }).is_err()
            {
                break;
            }
            current_qname = Some(qname.to_vec());
            group.push(ReadRec { id: read_id, record: record.clone(), segs: exons, decoded_seq });
        }
    }
    if !group.is_empty() {
        stats.read_groups += 1;
        batch.push((group_idx, group));
    }
    if !batch.is_empty() {
        let _ = tx.send(WorkItem { groups: batch });
    }
    // `tx` drops here, closing the channel.
    stats
}

struct WorkItem {
    groups: Vec<(usize, Vec<ReadRec>)>,
}

struct ResultItem {
    results: Vec<(usize, Result<Vec<hts_bam::Record>>)>,
}

fn process_group_records(
    group: &[ReadRec],
    g2t: &G2TTree,
    evaluator: &ReadEvaluator,
    fr: bool,
    rf: bool,
    paired_end: bool,
    ctx: &mut EvalContext,
) -> Result<Vec<hts_bam::Record>> {
    let mut output: Vec<hts_bam::Record> = Vec::with_capacity(group.len() * 2);

    let mut read_evals: Vec<ReadEval> = Vec::with_capacity(group.len());
    let shared_seq: Option<&[u8]> = group
        .iter()
        .find(|rec| !rec.record.is_unmapped() && !rec.decoded_seq.is_empty())
        .map(|rec| rec.decoded_seq.as_slice());

    for (idx, rec) in group.iter().enumerate() {
        if rec.record.is_unmapped() {
            continue;
        }

        let refid = if rec.record.tid() >= 0 { Some(rec.record.tid() as RefId) } else { None };
        if refid.is_none() {
            continue;
        }

        let mate_refid = if rec.record.mtid() >= 0 { Some(rec.record.mtid() as RefId) } else { None };

        let (strand, strand_from_tag) = splice_strand(&rec.record, fr, rf, paired_end);

        let cigar_ops: Vec<(u32, CigarKind)> = rec.record.cigar().iter()
            .map(hts_cigar_to_kind)
            .collect();
        let sequence: Option<&[u8]> = if rec.decoded_seq.is_empty() { None } else { Some(&rec.decoded_seq) };
        let name = rec.record.qname().to_vec();

        let read = ReadAln {
            strand,
            strand_from_tag,
            refid: refid.unwrap(),
            nh: 0,
            segs: rec.segs.clone(),
            cigar_ops,
            sequence: sequence.map(|s| s.to_vec()),
            name,
        };

        evaluator.evaluate(&read, rec.id, g2t, shared_seq, ctx);
        let matches: HashMap<Tid, ExonChainMatch> = ctx.matches
            .drain()
            .filter(|(_, m)| m.align.cigar.is_some())
            .collect();
        let read_len = rec.record.seq_len();
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

    // Count total output records (both read1 and read2 for paired), matching
    // C++ total_matches which increments once per r_align and once per m_align.
    let new_nh: u32 = output_groups.iter().map(|g| g.len() as u32).sum();
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
        let rec = &group[entry.record_idx];
        let out = build_projected_record(&rec.record, &rec.decoded_seq, &entry, nh, mapq, evaluator.long_reads)?;
        output.push(out);
    }

    Ok(output)
}

pub(crate) fn get_mapq(nh: u32, long_reads: bool) -> u8 {
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

pub(crate) fn assign_hit_indices(entries: &mut [OutputEntry]) {
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
}

fn get_alignment_start(record: &hts_bam::Record) -> Option<u32> {
    let pos = record.pos();
    if pos < 0 { None } else { Some((pos + 1) as u32) }
}

fn splice_strand(record: &hts_bam::Record, fr: bool, rf: bool, paired_end: bool) -> (char, bool) {
    if let Some(c) = get_char_tag(record, b"XS") {
        let strand = c as char;
        if strand == '+' || strand == '-' {
            return (strand, true);
        }
    }

    if let Some(c) = get_char_tag(record, b"ts") {
        let mut strand = c as char;
        if strand == '+' || strand == '-' {
            if record.is_reverse() {
                strand = if strand == '+' { '-' } else { '+' };
            }
            return (strand, true);
        }
    }
    // Match C++ behavior: if no splice-strand tags and no strand library flags,
    // return '.' (unspecified) so the evaluator checks both strands.
    // Only infer a definite strand when FR / RF is set.
    let rev = record.is_reverse();
    if fr || rf {
        // Treat read as paired if it has the SEGMENTED flag OR if --paired-end
        // was passed (C++: is_paired = brec->isPaired() || PAIRED_END).
        let is_paired = record.is_paired() || paired_end;
        let strand = if is_paired {
            if record.is_first_in_template() {
                if (rf && rev) || (fr && !rev) { '+' } else { '-' }
            } else if (rf && rev) || (fr && !rev) {
                '-'
            } else {
                '+'
            }
        } else if (rf && rev) || (fr && !rev) {
            '+'
        } else {
            '-'
        };
        (strand, false)
    } else {
        ('.', false)
    }
}

fn get_char_tag(record: &hts_bam::Record, tag: &[u8; 2]) -> Option<u8> {
    match record.aux(tag).ok()? {
        Aux::Char(c) => Some(c),
        Aux::String(s) => s.as_bytes().first().copied(),
        _ => None,
    }
}

fn get_mate_start(record: &hts_bam::Record) -> Option<u32> {
    if record.mtid() < 0 { return None; }
    let mpos = record.mpos();
    if mpos < 0 { None } else { Some((mpos + 1) as u32) }
}

fn get_hit_index(record: &hts_bam::Record) -> i32 {
    hts_aux_as_int(record.aux(b"HI").ok()).unwrap_or(0) as i32
}

pub(crate) fn find_mate_pairs(reads: &[ReadEval]) -> Vec<(usize, usize)> {
    let mut pending: HashMap<(u32, i32), usize> = HashMap::new();
    let mut pairs = Vec::new();

    for (idx, read) in reads.iter().enumerate() {
        if (read.flags & FLAG_PAIRED) == 0 || (read.flags & FLAG_MATE_UNMAPPED) != 0 {
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

pub(crate) fn assign_pair_order(i: usize, j: usize, reads: &[ReadEval]) -> (usize, usize) {
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

pub(crate) fn build_paired_groups(read: &ReadEval, mate: &ReadEval) -> Vec<Vec<OutputEntry>> {
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

pub(crate) fn build_unpaired_groups(read: &ReadEval) -> Vec<Vec<OutputEntry>> {
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

pub(crate) fn align_pos(match_info: &ExonChainMatch) -> u32 {
    if match_info.align.strand == '+' {
        match_info.align.fwpos
    } else {
        match_info.align.rcpos
    }
}

pub(crate) fn compute_template_length(
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
    record: &hts_bam::Record,
    decoded_seq: &[u8],
    entry: &OutputEntry,
    nh: u32,
    mapq: u8,
    long_reads: bool,
) -> Result<hts_bam::Record> {
    let mut out = record.clone();

    let original_reverse = record.is_reverse();
    let align_reverse = entry.align.align.is_reverse;
    let has_sequence = record.seq_len() > 0;
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

    // Build flags via u16 bit manipulation
    let mut flags = record.flags();
    // Clear: UNMAPPED, MATE_UNMAPPED, MATE_REVERSE, PROPER_PAIR, PAIRED, READ1, READ2
    flags &= !(FLAG_UNMAPPED | FLAG_MATE_UNMAPPED | FLAG_MATE_REVERSE | FLAG_PROPER_PAIR
                | FLAG_PAIRED | FLAG_READ1 | FLAG_READ2);

    if entry.align.align.primary_alignment {
        flags &= !FLAG_SECONDARY;
    } else {
        flags |= FLAG_SECONDARY;
    }

    if output_reverse {
        flags |= FLAG_REVERSE;
    } else {
        flags &= !FLAG_REVERSE;
    }

    let pos0 = if entry.align.align.strand == '+' {
        entry.align.align.fwpos
    } else {
        entry.align.align.rcpos
    };

    // Compute CIGAR
    let (mut cigar_ops, nm) = update_cigar_hts(record, entry.align.align.cigar.as_ref().unwrap())?;
    if reverse_action {
        cigar_ops.0.reverse();
    }

    // If reverse_action, we need to set seq+qual+cigar together
    if reverse_action {
        let mut seq = decoded_seq.to_vec();
        let raw_qual = record.qual();
        let qual_missing = raw_qual.first().copied() == Some(255);
        let mut qual: Vec<u8> = if qual_missing {
            vec![255u8; seq.len()]
        } else {
            raw_qual.to_vec()
        };
        reverse_complement(&mut seq);
        qual.reverse();
        out.set(record.qname(), Some(&cigar_ops), &seq, &qual);
    } else {
        out.set_cigar(Some(&cigar_ops));
    }

    out.set_tid(entry.tid as i32);
    out.set_pos(pos0 as i64);
    out.set_mapq(mapq);
    out.set_flags(flags);

    if let Some(mate) = entry.mate.as_ref() {
        flags |= FLAG_PAIRED | FLAG_PROPER_PAIR;
        if entry.is_first { flags |= FLAG_READ1; }
        if entry.is_last { flags |= FLAG_READ2; }
        if mate.is_reverse {
            flags |= FLAG_MATE_REVERSE;
        } else {
            flags &= !FLAG_MATE_REVERSE;
        }
        out.set_flags(flags);

        out.set_mtid(mate.tid as i32);
        out.set_mpos(mate.pos as i64);
        out.set_insert_size(compute_template_length(
            pos0,
            mate.pos,
            entry.read_len,
            mate.read_len,
            entry.is_first,
            entry.same_transcript,
        ) as i64);
    } else {
        flags &= !(FLAG_PAIRED | FLAG_PROPER_PAIR | FLAG_MATE_UNMAPPED
                    | FLAG_MATE_REVERSE | FLAG_READ1 | FLAG_READ2);
        out.set_flags(flags);
        out.set_mtid(-1);
        out.set_mpos(-1);
        out.set_insert_size(0);
    }

    // Extract AS from the *original* record before we modify aux on the clone.
    let gn_as = extract_alignment_score(record).unwrap_or(0) as f64;

    // Single-pass filter: remove SKIP_TAGS, preserve everything else.
    filter_aux_tags(&mut out);

    // Append our projected tags.  filter_aux_tags guarantees none of these
    // exist, so push_aux_unchecked is safe (skips the duplicate scan).
    out.push_aux_unchecked(b"NH", Aux::I32(nh as i32))?;
    out.push_aux_unchecked(b"HI", Aux::I32(entry.align.align.hit_index))?;
    out.push_aux_unchecked(b"XS", Aux::Char(entry.align.align.strand as u8))?;
    let score = if !long_reads {
        (1.0 + entry.align.align.similarity_score).powi(3) * 100.0
    } else {
        (gn_as + (entry.align.align.clip_score as f64)) * entry.align.align.similarity_score
    };
    out.push_aux_unchecked(b"AS", Aux::I32(score as i32))?;
    out.push_aux_unchecked(b"NM", Aux::I32(nm))?;

    Ok(out)
}


/// Single-pass aux tag filter: iterate raw aux bytes once, copy everything
/// except tags in `SKIP_TAGS`, then truncate the record's aux block and write
/// the filtered bytes back.  O(n) total instead of O(|SKIP_TAGS| * n).
fn filter_aux_tags(out: &mut hts_bam::Record) {
    let aux_start = out.inner.core.l_qname as usize
        + out.cigar_len() * std::mem::size_of::<u32>()
        + out.seq_len().div_ceil(2)
        + out.seq_len();

    let l_data = out.inner.l_data as usize;
    if aux_start >= l_data {
        return;
    }

    // Copy aux bytes so we can read while mutating the record.
    let aux_bytes: Vec<u8> = unsafe {
        std::slice::from_raw_parts(out.inner.data.add(aux_start), l_data - aux_start)
    }
    .to_vec();

    let mut kept = Vec::with_capacity(aux_bytes.len());
    let mut pos = 0;
    while pos + 3 <= aux_bytes.len() {
        let tag = [aux_bytes[pos], aux_bytes[pos + 1]];
        let type_byte = aux_bytes[pos + 2];
        let val_len = match type_byte {
            b'A' | b'c' | b'C' => 1,
            b's' | b'S' => 2,
            b'i' | b'I' | b'f' => 4,
            b'd' => 8,
            b'Z' | b'H' => match aux_bytes[pos + 3..].iter().position(|&b| b == 0) {
                Some(end) => end + 1,
                None => break,
            },
            b'B' => {
                if pos + 8 > aux_bytes.len() {
                    break;
                }
                let subtype = aux_bytes[pos + 3];
                let count = u32::from_le_bytes([
                    aux_bytes[pos + 4],
                    aux_bytes[pos + 5],
                    aux_bytes[pos + 6],
                    aux_bytes[pos + 7],
                ]) as usize;
                let elem_size = match subtype {
                    b'c' | b'C' => 1,
                    b's' | b'S' => 2,
                    b'i' | b'I' | b'f' => 4,
                    b'd' => 8,
                    _ => break,
                };
                5 + count * elem_size
            }
            _ => break,
        };
        let entry_len = 3 + val_len; // tag(2) + type(1) + value
        if pos + entry_len > aux_bytes.len() {
            break;
        }

        let should_skip = SKIP_TAGS.iter().any(|s| *s == tag);
        if !should_skip {
            kept.extend_from_slice(&aux_bytes[pos..pos + entry_len]);
        }
        pos += entry_len;
    }

    // Truncate aux, then write back kept bytes.
    // kept.len() <= original aux size, so it fits within the already-allocated m_data.
    out.inner.l_data = aux_start as i32;
    if !kept.is_empty() {
        unsafe {
            std::ptr::copy_nonoverlapping(
                kept.as_ptr(),
                out.inner.data.add(aux_start),
                kept.len(),
            );
        }
        out.inner.l_data = (aux_start + kept.len()) as i32;
    }
}

fn extract_alignment_score(record: &hts_bam::Record) -> Option<i32> {
    hts_aux_as_int(record.aux(b"AS").ok()).map(|v| v as i32)
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


/// Extract an integer value from any numeric Aux variant.
fn hts_aux_as_int(aux: Option<Aux<'_>>) -> Option<i64> {
    match aux? {
        Aux::I8(v)  => Some(v as i64),
        Aux::U8(v)  => Some(v as i64),
        Aux::I16(v) => Some(v as i64),
        Aux::U16(v) => Some(v as i64),
        Aux::I32(v) => Some(v as i64),
        Aux::U32(v) => Some(v as i64),
        _ => None,
    }
}

/// Convert a rust-htslib Cigar operation to a `(length, CigarKind)` pair.
fn hts_cigar_to_kind(op: &hts_bam::record::Cigar) -> (u32, CigarKind) {
    use hts_bam::record::Cigar::*;
    match op {
        Match(n)    => (*n, CigarKind::Match),
        Ins(n)      => (*n, CigarKind::Insertion),
        Del(n)      => (*n, CigarKind::Deletion),
        RefSkip(n)  => (*n, CigarKind::Skip),
        SoftClip(n) => (*n, CigarKind::SoftClip),
        HardClip(n) => (*n, CigarKind::HardClip),
        Pad(n)      => (*n, CigarKind::Pad),
        Equal(n)    => (*n, CigarKind::SequenceMatch),
        Diff(n)     => (*n, CigarKind::SequenceMismatch),
    }
}

/// Convert a rust-htslib Cigar operation to a noodles SamCigarOp.
fn hts_cigar_to_sam_op(op: &hts_bam::record::Cigar) -> SamCigarOp {
    let (len, kind) = hts_cigar_to_kind(op);
    SamCigarOp::new(kind, len as usize)
}

fn update_cigar_hts(
    record: &hts_bam::Record,
    ideal: &crate::evaluate::Cigar,
) -> Result<(CigarString, i32)> {
    let real_ops: Vec<SamCigarOp> = record.cigar().iter()
        .map(hts_cigar_to_sam_op)
        .collect();
    let (sam_cigar, nm) = update_cigar_ops(&real_ops, ideal);
    // Convert SamCigar to htslib CigarString
    let hts_ops: Vec<hts_bam::record::Cigar> = sam_cigar.as_ref().iter()
        .map(|op| {
            let kind = op.kind();
            let len = op.len();
            use hts_bam::record::Cigar::*;
            let n = len as u32;
            match kind {
                CigarKind::Match => Match(n),
                CigarKind::Insertion => Ins(n),
                CigarKind::Deletion => Del(n),
                CigarKind::Skip => RefSkip(n),
                CigarKind::SoftClip => SoftClip(n),
                CigarKind::HardClip => HardClip(n),
                CigarKind::Pad => Pad(n),
                CigarKind::SequenceMatch => Equal(n),
                CigarKind::SequenceMismatch => Diff(n),
            }
        })
        .collect();
    Ok((CigarString(hts_ops), nm))
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

// ── Library API helpers ────────────────────────────────────────────────────────

/// Convert a SAM op-code byte (0–8) to a noodles `CigarKind`.
pub(crate) fn sam_op_to_kind(op: u8) -> Option<CigarKind> {
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
pub(crate) fn segs_from_ops(ref_start: i64, ops: &[(u32, CigarKind)]) -> Vec<alignment::Segment> {
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


