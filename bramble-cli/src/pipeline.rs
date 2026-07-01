// pipeline.rs: the bramble-cli BAM read → project → BAM write driver. The pure
// projection logic lives in the `bramble-rs` library; this file owns the htslib
// I/O, threading, and SAM-record construction around it.
#![allow(dead_code)]
use crate::alignment;
use crate::bam_input::BamInput;
use crate::cli::Args;
use bramble_rs::evaluate::{EvalContext, ExonChainMatch, ReadAln, ReadEvaluator};
use bramble_rs::g2t::G2TTree;
use bramble_rs::groups::{
    ReadEval, OutputEntry, find_mate_pairs, assign_pair_order, build_paired_groups,
    build_unpaired_groups, assign_hit_indices, compute_template_length,
    FLAG_PAIRED, FLAG_PROPER_PAIR, FLAG_UNMAPPED, FLAG_MATE_UNMAPPED, FLAG_REVERSE,
    FLAG_MATE_REVERSE, FLAG_READ1, FLAG_READ2, FLAG_SECONDARY,
};
use bramble_rs::types::{HashMap, ReadId, RefId, Tid};
use anyhow::Result;
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use noodles::sam::alignment::record::cigar::{op::Kind as CigarKind, Op as SamCigarOp};
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
/// or will be replaced with bramble-computed values).
const SKIP_TAGS: &[[u8; 2]] = &[
    *b"NH", *b"HI", *b"AS", *b"CG",
    *b"XS", *b"SA",
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


pub fn run(
    args: &Args,
    hts_header: &rust_htslib::bam::Header,
    bam: &mut BamInput,
    g2t: &G2TTree,
    fasta: Option<&bramble_rs::fasta::FastaDb>,
) -> Result<Stats> {
    // BGZF (de)compression of the input and (large) output BAMs is the dominant
    // cost. Share ONE htslib thread pool between the reader and the writer (via
    // hts_set_thread_pool) so input decompression and output compression are both
    // parallelized without spinning up two separate per-file pools. Declared
    // before `writer` so it outlives it (the writer flushes on drop using it).
    let io_threads = (args.threads.get() as usize).max(4);
    let tpool = rust_htslib::tpool::ThreadPool::new(io_threads as u32)?;
    let mut writer = hts_bam::Writer::from_path(&args.out_bam, hts_header, hts_bam::Format::Bam)?;
    writer.set_thread_pool(&tpool)?;
    bam.reader.set_thread_pool(&tpool)?;

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

    let evaluator = ReadEvaluator {
        lr:                   args.lr,
        lr_hq:                args.lr_hq,
        strict:               args.strict,
        use_fasta:            fasta.is_some(),
        max_clip:             args.max_clip,
        max_ins:              args.max_ins,
        max_junc_gap:         args.max_junc_gap,
        similarity_threshold: args.similarity_threshold,
        small_exon_size:      args.small_exon_size,
        junc_miss_discount:   None,
    };
    
    let fr = args.fr;
    let rf = args.rf;
    let worker_count = args.threads.get() as usize;
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
                                group, g2t_ref, evaluator_ref, fr, rf, &mut ctx,
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
                        group, g2t_ref, evaluator_ref, fr, rf, &mut ctx,
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
                                group, g2t_ref, evaluator_ref, fr, rf, &mut ctx,
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

        let (strand, strand_from_tag) = splice_strand(&rec.record, fr, rf);

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
        let (r_idx, m_idx) = assign_pair_order(i, j, &read_evals);
        let read = &read_evals[r_idx];
        let mate = &read_evals[m_idx];

        // Match C++ process_mate_pair logic. In C++, the lower-indexed record
        // (first in BAM order, typically R1) is `this_read` and its mate is
        // `mate_read`. If this_read is null (no matches), the function returns
        // immediately and the mate's matches are lost. If mate_read is null,
        // this_read is emitted as unpaired.
        //
        // assign_pair_order makes R1 = "read" and R2 = "mate", matching C++.
        if read.matches.is_empty() {
            // R1 has no matches → C++ returns immediately, both dropped
        } else if mate.matches.is_empty() {
            // R2 has no matches → C++ emits R1 as unpaired
            output_groups.extend(build_unpaired_groups(read));
        } else {
            // Both have matches → try pairing; if it fails (case 5: no common
            // transcripts, at least one multi-mapped), both are dropped
            output_groups.extend(build_paired_groups(read, mate));
        }
        paired[i] = true;
        paired[j] = true;
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
        let nh = new_nh;
        let mapq = get_mapq(nh, evaluator.lr);
        let rec = &group[entry.record_idx];
        let out = build_projected_record(&rec.record, &rec.decoded_seq, &entry, nh, mapq, evaluator.lr)?;
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

fn get_alignment_start(record: &hts_bam::Record) -> Option<u32> {
    let pos = record.pos();
    if pos < 0 { None } else { Some((pos + 1) as u32) }
}

fn splice_strand(record: &hts_bam::Record, fr: bool, rf: bool) -> (char, bool) {
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
        let is_paired = record.is_paired();
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
    let has_sequence = record.seq_len() > 0;
    // C++ approach: reverse_complement_bam is called when strand == '-',
    // which XORs FLAG_REVERSE and reverses seq/qual/cigar (if sequence present).
    // This applies identically for both short and long reads.
    let reverse_action = entry.align.align.strand == '-';
    let output_reverse = original_reverse ^ reverse_action;

    // Build flags via u16 bit manipulation
    let mut flags = record.flags();
    // Original mate identity (PAIRED/READ1/READ2), preserved for orphans below so
    // the projected record keeps its read1/read2 designation like C++ (which only
    // adjusts SECONDARY + the mate fields and never strips mate identity).
    let orig_mate_id = flags & (FLAG_PAIRED | FLAG_READ1 | FLAG_READ2);
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
    let (mut cigar_ops, _nm) = update_cigar_hts(record, entry.align.align.cigar.as_ref().unwrap())?;
    if reverse_action {
        cigar_ops.0.reverse();
    }

    // If reverse_action, we need to set seq+qual+cigar together.
    // If sequence is absent (l_qseq <= 0), C++ just XORs the flag without
    // touching seq/qual — we only need to set the cigar.
    if reverse_action && has_sequence {
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
        flags |= FLAG_PAIRED;
        if entry.same_transcript {
            flags |= FLAG_PROPER_PAIR;
        } else {
            flags &= !FLAG_PROPER_PAIR;
        }
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
            entry.same_transcript,
        ) as i64);
    } else {
        // Orphan: this read's mate did not project onto this transcript. Preserve
        // the original PAIRED/READ1/READ2 designation (C++ keeps mate identity),
        // but mark the mate unmapped and clear the (genomic) mate coordinates,
        // since there is no mate placement in transcriptome space.
        flags &= !(FLAG_PROPER_PAIR | FLAG_MATE_REVERSE);
        flags |= orig_mate_id;
        if flags & FLAG_PAIRED != 0 {
            flags |= FLAG_MATE_UNMAPPED;
        }
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
    // Projected alignment score, matching C++ exactly. C++ (core.cpp) only calls
    // set_as_tag — `(genome_AS + clip_score) * similarity_score` — for LONG reads;
    // for SHORT reads it leaves the original genome (STAR) AS untouched. So short
    // reads pass the genome AS through verbatim (per-read uniform), and only long
    // reads fold in the projection similarity. The previous Rust behavior recomputed
    // AS for short reads too — first as `(1+similarity)^3*100`, then as the long-read
    // formula — either of which applies the junc_hits-driven similarity variance that
    // C++ never applies to short reads, so downstream quantifiers split multi-isoform
    // reads unevenly on score noise (dropping accuracy vs C++).
    let score = if long_reads {
        (gn_as + (entry.align.align.clip_score as f64)) * entry.align.align.similarity_score
    } else {
        gn_as
    };
    out.push_aux_unchecked(b"AS", Aux::I32(score as i32))?;

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

        let should_skip = SKIP_TAGS.contains(&tag);
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
    ideal: &bramble_rs::evaluate::Cigar,
) -> Result<(CigarString, i32)> {
    let real_ops: Vec<SamCigarOp> = record.cigar().iter()
        .map(hts_cigar_to_sam_op)
        .collect();
    let (sam_cigar, nm) = bramble_rs::cigar::update_cigar_ops(&real_ops, ideal);
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
