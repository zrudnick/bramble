use crate::annotation::{Exon, Transcript};
use crate::evaluate::{ExonStatus, ReadEvaluationConfig};
use crate::fasta::FastaDb;
use crate::types::{RefId, Tid};
use anyhow::{anyhow, Result};
use coitrees::{BasicCOITree, Interval, IntervalTree as CoitreeIntervalTree};
use std::collections::{HashMap, HashSet};

#[derive(Debug, Clone)]
pub struct IntervalData {
    pub start: u32,
    pub end: u32,
    pub idx: u8,
    pub pos_start: u32,
    pub seq: String,

    pub has_prev: bool,
    pub has_next: bool,
    pub prev_start: u32,
    pub prev_end: u32,
    pub next_start: u32,
    pub next_end: u32,

    pub transcript_len: u32,
}

#[derive(Debug, Clone)]
pub struct GuideExon {
    pub start: u32,
    pub end: u32,
    pub pos: u32,
    pub pos_start: u32,
    pub exon_id: u8,
    pub seq: String,
    pub left_ins: i32,
    pub right_ins: i32,
    pub left_gap: i32,
    pub right_gap: i32,

    pub has_prev: bool,
    pub has_next: bool,
    pub prev_start: u32,
    pub prev_end: u32,
    pub next_start: u32,
    pub next_end: u32,

    pub transcript_len: u32,
}

#[derive(Debug, Clone, Default)]
pub struct IITData {
    pub tid: Tid,
    pub exon_id: u8,
    pub pos_start: u32,
    pub has_prev: bool,
    pub has_next: bool,
    pub prev_start: u32,
    pub prev_end: u32,
    pub next_start: u32,
    pub next_end: u32,
    pub transcript_len: u32,
}

pub struct IntervalTree {
    intervals: Vec<Interval<IITData>>,
    tree: Option<BasicCOITree<IITData, u32>>,
    seqs: HashMap<(Tid, u8), String>,
}

impl IntervalTree {
    pub fn new() -> Self {
        Self {
            intervals: Vec::new(),
            tree: None,
            seqs: HashMap::new(),
        }
    }

    pub fn add_interval(&mut self, tid: Tid, interval: &IntervalData) {
        let data = IITData {
            tid,
            exon_id: interval.idx,
            pos_start: interval.pos_start,
            has_prev: interval.has_prev,
            has_next: interval.has_next,
            prev_start: interval.prev_start,
            prev_end: interval.prev_end,
            next_start: interval.next_start,
            next_end: interval.next_end,
            transcript_len: interval.transcript_len,
        };
        if !interval.seq.is_empty() {
            self.seqs.insert((tid, interval.idx), interval.seq.clone());
        }
        // COITree intervals are end-inclusive; convert [start, end) -> [start, end-1].
        let first = interval.start as i32;
        let last = interval.end.saturating_sub(1) as i32;
        if last >= first {
            self.intervals.push(Interval::new(first, last, data));
        }
    }

    pub fn index(&mut self) {
        self.tree = Some(BasicCOITree::new(&self.intervals));
    }

    pub fn find_overlapping_for_tid(
        &self,
        qstart: u32,
        qend: u32,
        tid: Tid,
    ) -> Option<GuideExon> {
        let tree = self.tree.as_ref()?;
        if qstart == 0 && qend == 0 {
            return None;
        }

        let q_last = qend.saturating_sub(1) as i32;
        let mut found: Option<GuideExon> = None;

        tree.query(qstart as i32, q_last, |node| {
            if found.is_some() {
                return;
            }
            if node.metadata.tid == tid {
                let s = node.first.max(0) as u32;
                let e_excl = (node.last + 1).max(0) as u32;
                let data = &node.metadata;
                found = Some(GuideExon {
                    start: s,
                    end: e_excl,
                    pos: 0,
                    pos_start: data.pos_start,
                    exon_id: data.exon_id,
                    seq: self
                        .seqs
                        .get(&(data.tid, data.exon_id))
                        .cloned()
                        .unwrap_or_default(),
                    left_ins: 0,
                    right_ins: 0,
                    left_gap: 0,
                    right_gap: 0,
                    has_prev: data.has_prev,
                    has_next: data.has_next,
                    prev_start: data.prev_start,
                    prev_end: data.prev_end,
                    next_start: data.next_start,
                    next_end: data.next_end,
                    transcript_len: data.transcript_len,
                });
            }
        });

        found
    }

    pub fn find_overlapping(
        &self,
        qstart: u32,
        qend: u32,
        strand: char,
        config: &ReadEvaluationConfig,
        status: ExonStatus,
    ) -> HashMap<Tid, GuideExon> {
        let mut exons: HashMap<Tid, GuideExon> = HashMap::new();
        let tree = match self.tree.as_ref() {
            Some(t) => t,
            None => return exons,
        };

        let q_last = qend.saturating_sub(1) as i32;

        tree.query(qstart as i32, q_last, |node| {
            let s = node.first.max(0) as u32;
            let e_excl = (node.last + 1).max(0) as u32;

            if (!config.soft_clips || config.strict) && (qstart < s || e_excl < qend) {
                return;
            }

            let data = &node.metadata;
            let mut exon = GuideExon {
                start: s,
                end: e_excl,
                pos: 0,
                pos_start: data.pos_start,
                exon_id: data.exon_id,
                seq: self
                    .seqs
                    .get(&(data.tid, data.exon_id))
                    .cloned()
                    .unwrap_or_default(),
                left_ins: 0,
                right_ins: 0,
                left_gap: 0,
                right_gap: 0,
                has_prev: data.has_prev,
                has_next: data.has_next,
                prev_start: data.prev_start,
                prev_end: data.prev_end,
                next_start: data.next_start,
                next_end: data.next_end,
                transcript_len: data.transcript_len,
            };

            let mut rejected = false;

            if strand == '+' {
                if s <= qstart {
                    exon.pos = (qstart - s) + data.pos_start;
                    exon.left_gap = (qstart - s) as i32;
                    if matches!(status, ExonStatus::MiddleExon | ExonStatus::LastExon)
                        && (exon.left_gap as u32) > config.max_junc_gap
                    {
                        rejected = true;
                    }
                } else {
                    exon.pos = data.pos_start;
                    exon.left_ins = (s - qstart) as i32;
                    if matches!(status, ExonStatus::MiddleExon | ExonStatus::LastExon) {
                        if (exon.left_ins as u32) > config.max_ins {
                            rejected = true;
                        }
                    } else if (exon.left_ins as u32) > config.max_clip {
                        rejected = true;
                    }
                }

                if !rejected {
                    if e_excl < qend {
                        exon.right_ins = (qend - e_excl) as i32;
                        if matches!(status, ExonStatus::FirstExon | ExonStatus::MiddleExon) {
                            if (exon.right_ins as u32) > config.max_ins {
                                rejected = true;
                            }
                        } else if (exon.right_ins as u32) > config.max_clip {
                            rejected = true;
                        }
                    } else if qend < e_excl {
                        exon.right_gap = (e_excl - qend) as i32;
                        if matches!(status, ExonStatus::FirstExon | ExonStatus::MiddleExon)
                            && (exon.right_gap as u32) > config.max_junc_gap
                        {
                            rejected = true;
                        }
                    }
                }
            } else {
                if qend <= e_excl {
                    exon.pos = (e_excl - qend) + data.pos_start;
                    exon.right_gap = (e_excl - qend) as i32;
                    if matches!(status, ExonStatus::FirstExon | ExonStatus::MiddleExon)
                        && (exon.right_gap as u32) > config.max_junc_gap
                    {
                        rejected = true;
                    }
                } else {
                    exon.pos = data.pos_start;
                    exon.right_ins = (qend - e_excl) as i32;
                    // NOTE: Match C++ behavior (g2t.cpp) which effectively
                    // treats this condition as always true due to `|| MIDDLE_EXON`.
                    // That means `max_ins` is enforced for all statuses here.
                    if (exon.right_ins as u32) > config.max_ins {
                        rejected = true;
                    }
                }

                if !rejected {
                    if qstart < s {
                        exon.left_ins = (s - qstart) as i32;
                        if matches!(status, ExonStatus::MiddleExon | ExonStatus::LastExon) {
                            if (exon.left_ins as u32) > config.max_ins {
                                rejected = true;
                            }
                        } else if (exon.left_ins as u32) > config.max_clip {
                            rejected = true;
                        }
                    } else if s < qstart {
                        exon.left_gap = (qstart - s) as i32;
                        if matches!(status, ExonStatus::MiddleExon | ExonStatus::LastExon)
                            && (exon.left_gap as u32) > config.max_junc_gap
                        {
                            rejected = true;
                        }
                    }
                }
            }

            if !rejected {
                exons.insert(data.tid, exon);
            }
        });

        exons
    }
}

impl Default for IntervalTree {
    fn default() -> Self {
        Self::new()
    }
}

pub struct G2TTree {
    // Per refid: (forward, reverse)
    trees: HashMap<RefId, (IntervalTree, IntervalTree)>,
    name_id_map: HashMap<Tid, String>,
    tid_names: Vec<String>,
}

impl G2TTree {
    pub fn new() -> Self {
        Self {
            trees: HashMap::new(),
            name_id_map: HashMap::new(),
            tid_names: Vec::new(),
        }
    }

    fn get_trees_mut(&mut self, refid: RefId) -> &mut (IntervalTree, IntervalTree) {
        self.trees.entry(refid).or_insert_with(|| {
            (IntervalTree::new(), IntervalTree::new())
        })
    }

    fn get_tree_mut(&mut self, refid: RefId, strand: char) -> &mut IntervalTree {
        let (fw, rc) = self.get_trees_mut(refid);
        if strand == '-' {
            rc
        } else {
            fw
        }
    }

    pub fn insert_tid_name(&mut self, tid: Tid, tid_name: &str) {
        if let std::collections::hash_map::Entry::Vacant(entry) = self.name_id_map.entry(tid) {
            entry.insert(tid_name.to_string());
            self.tid_names.push(tid_name.to_string());
        }
    }

    pub fn add_interval(
        &mut self,
        refid: RefId,
        tid: Tid,
        interval: &IntervalData,
        strand: char,
    ) {
        let tree = self.get_tree_mut(refid, strand);
        tree.add_interval(tid, interval);
    }

    pub fn index_ref(&mut self, refid: RefId) {
        let (fw, rc) = self.get_trees_mut(refid);
        fw.index();
        rc.index();
    }

    pub fn get_guide_exon_for_tid(
        &self,
        refid: RefId,
        strand: char,
        tid: Tid,
        start: u32,
        end: u32,
    ) -> Option<GuideExon> {
        let tree = self.get_tree(refid, strand)?;
        tree.find_overlapping_for_tid(start, end, tid)
    }

    pub fn get_guide_exons(
        &self,
        refid: RefId,
        strand: char,
        start: u32,
        end: u32,
        config: &ReadEvaluationConfig,
        status: ExonStatus,
    ) -> HashMap<Tid, GuideExon> {
        if let Some(tree) = self.get_tree(refid, strand) {
            tree.find_overlapping(start, end, strand, config, status)
        } else {
            HashMap::new()
        }
    }

    fn get_tree(&self, refid: RefId, strand: char) -> Option<&IntervalTree> {
        let pair = self.trees.get(&refid)?;
        if strand == '-' {
            Some(&pair.1)
        } else {
            Some(&pair.0)
        }
    }
}

impl Default for G2TTree {
    fn default() -> Self {
        Self::new()
    }
}

pub fn build_g2t(
    transcripts: &[Transcript],
    refname_to_id: &HashMap<String, RefId>,
    fasta: Option<&FastaDb>,
) -> Result<G2TTree> {
    let mut g2t = G2TTree::new();
    let mut touched_refids: HashSet<RefId> = HashSet::new();

    // Assign tid incrementally in stable order of transcripts slice
    for (tid_idx, tx) in transcripts.iter().enumerate() {
        let tid = tid_idx as Tid;
        g2t.insert_tid_name(tid, &tx.id);

        let refid = *refname_to_id
            .get(&tx.seqname)
            .ok_or_else(|| anyhow!("reference '{}' not found in BAM header", tx.seqname))?;
        touched_refids.insert(refid);

        let mut exons: Vec<Exon> = tx.exons.clone();
        exons.sort_by_key(|e| e.start);

        let mut pos_start: u32 = 0;
        let mut intervals: Vec<IntervalData> = Vec::with_capacity(exons.len());
        let exon_iter: Box<dyn Iterator<Item = (usize, &Exon)>> = if tx.strand == '-' {
            Box::new(exons.iter().enumerate().rev())
        } else {
            Box::new(exons.iter().enumerate())
        };

        for (idx, exon) in exon_iter {
            let len = exon.end.saturating_sub(exon.start);
            let seq = if let Some(fa) = fasta {
                fa.get_slice(&tx.seqname, exon.start, exon.end)
                    .map(|v| String::from_utf8_lossy(&v).to_string())
                    .unwrap_or_default()
            } else {
                String::new()
            };
            intervals.push(IntervalData {
                start: exon.start,
                end: exon.end,
                idx: idx as u8,
                pos_start,
                seq,
                has_prev: false,
                has_next: false,
                prev_start: 0,
                prev_end: 0,
                next_start: 0,
                next_end: 0,
                transcript_len: 0,
            });
            pos_start = pos_start.saturating_add(len);
        }

        let transcript_len = pos_start;

        for k in 0..intervals.len() {
            let mut interval = intervals[k].clone();
            if k > 0 {
                interval.has_prev = true;
                interval.prev_start = intervals[k - 1].start;
                interval.prev_end = intervals[k - 1].end;
            }
            if k + 1 < intervals.len() {
                interval.has_next = true;
                interval.next_start = intervals[k + 1].start;
                interval.next_end = intervals[k + 1].end;
            }
            interval.transcript_len = transcript_len;

            g2t.add_interval(refid, tid, &interval, tx.strand);
        }

    }

    // Index all trees only after all intervals are inserted.
    for refid in touched_refids {
        g2t.index_ref(refid);
    }

    Ok(g2t)
}
