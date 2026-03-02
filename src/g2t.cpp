
#include <numeric>
#include <algorithm>
#include <iostream>
#include <memory>
#include <vector>

#include <IITree.h>
#include "types.h"
#include "bramble.h"
#include "evaluate.h"
#include "g2t.h"

extern bool SOFT_CLIPS;
extern bool STRICT;
extern bool USE_FASTA;
extern bool LONG_READS;

extern GFastaDb *gfasta;  // FASTA file, -S

namespace bramble {

  IITData::IITData() {}

  IntervalTree::IntervalTree() {
    iit = new IITree<int, IITData>();
  }

  IntervalTree::~IntervalTree() { 
    delete(iit); 
  }

  // Add guide exon with TID and transcript start
  void IntervalTree::addInterval(const tid_t &tid, IntervalData interval, 
                                const char* ref_name) {
    IITData node = IITData();
    node.tid = tid;
    node.exon_id = interval.idx;
    node.pos_start = interval.pos_start;

    node.has_prev = interval.has_prev;
    node.has_next = interval.has_next;
    node.prev_start = interval.prev_start;
    node.prev_end = interval.prev_end;
    node.next_start = interval.next_start;
    node.next_end = interval.next_end;
    node.transcript_len = interval.transcript_len;

    // grab here to avoid multithreading issues
    if (USE_FASTA) {
      GFaSeqGet* faseq = gfasta->fetch(ref_name);
      const char *seq = faseq->copyRange(interval.start, 
        interval.end - 1, false, true);
      node.seq = seq;
    }

    iit->add(interval.start, interval.end, node);
  }

  void IntervalTree::indexTree() {
    iit->index();
  }

  struct iitData {
    uint32_t start;
    uint32_t end;
    IITData data;
  };

  bool
  IntervalTree::findOverlappingForTid(uint32_t qstart, uint32_t qend, tid_t tid,
                                      GuideExon &gexon) {
    
    if (qstart == 0 && qend == 0) return false;
    std::vector<size_t> hits;
    iit->overlap(qstart, qend, hits);

    for (size_t idx : hits) {
      uint32_t s = iit->start(idx);
      uint32_t e = iit->end(idx);

      const IITData& data = iit->data(idx);
      if (data.tid == tid) {
        gexon.start = s;
        gexon.end = e;
        gexon.transcript_len = data.transcript_len;
        gexon.pos_start = data.pos_start;
        gexon.exon_id = data.exon_id;
        gexon.has_prev = data.has_prev;
        gexon.has_next = data.has_next;
        gexon.prev_start = data.prev_start;
        gexon.prev_end = data.prev_end;
        gexon.next_start = data.next_start;
        gexon.next_end = data.next_end;
        gexon.seq = data.seq; // may have to be cleared
        return true;
      }
    }

    return false;
  }

  bool
  IntervalTree::findOverlapping(uint32_t qstart, uint32_t qend, 
                                char strand, ReadEvaluationConfig config, 
                                ExonStatus status, 
                                std::vector<GuideExon> &gexons) {
    std::vector<size_t> hits;
    
    iit->overlap(qstart, qend, hits);
    if (hits.empty()) return false;

    uint32_t pos;
    uint32_t left_gap;
    uint32_t left_ins;
    uint32_t right_gap;
    uint32_t right_ins;

    // remove tids that appear >1
    // means that junctions have been ignored

    // update: should be handled in create_match loop
    // if (LONG_READS) { // won't happen with strict short read parameters
    //   thread_local unordered_set<tid_t> seen;
    //   seen.clear();
    //   thread_local unordered_set<tid_t> duplicates;
    //   duplicates.clear();

    //   for (size_t idx : hits) {
    //     tid_t tid = iit->data(idx).tid;
    //     if (!seen.emplace(tid).second) {
    //       duplicates.emplace(tid);
    //     }
    //   }

    //   hits.erase(
    //     std::remove_if(hits.begin(), hits.end(), [&](size_t idx) {
    //       return duplicates.count(iit->data(idx).tid);
    //     }),
    //     hits.end()
    //   );
    // }

    for (size_t idx : hits) {
      uint32_t s = iit->start(idx);
      uint32_t e = iit->end(idx);

      pos = 0;
      left_gap = 0;
      left_ins = 0;
      right_gap = 0;
      right_ins = 0;

      const IITData& data = iit->data(idx);

      if (strand == '+') {
        if (s <= qstart) {
          // Left junction gap: query starts after interval starts
          pos = (qstart - s) + data.pos_start;
          left_gap = qstart - s;
          if (status == MIDDLE_EXON || status == LAST_EXON) {
            if (left_gap > config.max_junc_gap) continue;
          }
        } else {
          pos = data.pos_start;
          // Left clip: query extends before interval
          left_ins = s - qstart;
          if (status == MIDDLE_EXON || status == LAST_EXON) {
            if (left_ins > config.max_ins) continue;
          } else {
            if (left_ins > config.max_clip) continue;
          }
        }

        // Check right side
        if (e < qend) {
          // Right clip: query extends past interval
          right_ins = qend - e;
          if (status == FIRST_EXON || status == MIDDLE_EXON) {
            if (right_ins > config.max_ins) continue;
          } else {
            if (right_ins > config.max_clip) continue;
          }
        } else if (qend < e) {
          // Right junction gap: query ends before interval ends
          right_gap = e - qend;
          if (status == FIRST_EXON || status == MIDDLE_EXON) {
            if (right_gap > config.max_junc_gap) continue;
          }
        }

      } else {  // strand == '-'
        if (qend <= e) {
          pos = (e - qend) + data.pos_start;
          right_gap = e - qend;
          // Right junction gap: query ends before interval ends
          if (status == FIRST_EXON || status == MIDDLE_EXON) {
            if (right_gap > config.max_junc_gap) continue;
          }
        } else {
          pos = data.pos_start;
          // Right clip: query extends past interval
          right_ins = qend - e;
          if (status == FIRST_EXON || MIDDLE_EXON) {
            if (right_ins > config.max_ins) continue;
          } else {
            if (right_ins > config.max_clip) continue;
          }
        }

        // Check left side
        if (qstart < s) {
          // Left clip: query extends before interval
          left_ins = s - qstart;
          if (status == MIDDLE_EXON || status == LAST_EXON) {
            if (left_ins > config.max_ins) continue;
          } else {
            if (left_ins > config.max_clip) continue;
          }
        } else if (s < qstart) {
          left_gap = qstart - s;
          // Left junction gap: query starts after interval starts
          if (status == MIDDLE_EXON || status == LAST_EXON) {
            if (left_gap > config.max_junc_gap) continue;
          }
        }
      }

      auto exon = GuideExon();
      exon.start = s;
      exon.end = e;
      exon.tid = data.tid;

      exon.pos = pos;
      exon.left_gap = left_gap;
      exon.left_ins = left_ins;
      exon.right_gap = right_gap;
      exon.right_ins = right_ins;

      exon.pos_start = data.pos_start;
      exon.exon_id = data.exon_id;

      exon.has_prev = data.has_prev;
      exon.has_next = data.has_next;
      exon.prev_start = data.prev_start;
      exon.prev_end = data.prev_end;
      exon.next_start = data.next_start;
      exon.next_end = data.next_end;

      exon.transcript_len = data.transcript_len;
      exon.seq = data.seq; // may have to be cleared

      gexons.emplace_back(exon);
    }

    return (!gexons.empty());
  }

  g2tTree::g2tTree() = default;
    
  g2tTree::~g2tTree() {
    // Clean up all allocated trees
    for (auto& pair : trees) {
      delete pair.first;   // forward tree
      delete pair.second;  // reverse tree
    }
    trees.clear();
  }

  // Get or create tree pair for a given refid
  void g2tTree::createTree(int refid) {
    // Create new tree pair for this refid
    IntervalTree* fw_tree = new IntervalTree();
    IntervalTree* rc_tree = new IntervalTree();
    trees.emplace_back(std::make_pair(fw_tree, rc_tree));
  }

  IntervalTree *g2tTree::getTreeForStrand(int refid, char strand) {
    if (refid < 0 || (size_t)refid >= trees.size()) return nullptr;
    auto pair = trees[refid];
    if (strand == '+' || strand == 1)
      return pair.first;
    else if (strand == '-' || strand == -1)
      return pair.second;
    else
      return nullptr;
  }

  // If tid_name is already part of the tree, then return the existing
  // tid for it, otherwise give it the next available tid and return that.
  tid_t g2tTree::insertTidString(const char *&tid_name, BamIO* io) {
    tid_t tid = io->get_tid(tid_name);

    auto lookup = name_id_map.find(tid);
    if (lookup == name_id_map.end()) {
      name_id_map.insert({tid, tid_name});
      tid_names.push_back(tid_name);
    }

    return tid;
  }

  // get the string name given an id
  const std::string &g2tTree::getTidName(tid_t id) {
    return (name_id_map.find(id) != name_id_map.end()) ? 
      name_id_map[id] : invalid_name;
  }

  // Add guide exon with TID and transcript start
  void g2tTree::addInterval(refid_t refid, const tid_t &tid, 
                            IntervalData interval, char strand,
                            const char* ref_name) {
    IntervalTree *tree = getTreeForStrand(refid, strand);
    tree->addInterval(tid, interval, ref_name);
  }

  // Index tree after guides have been added
  void g2tTree::indexTrees(int refid) {
    IntervalTree *fw_tree = getTreeForStrand(refid, '+');
    IntervalTree *rc_tree = getTreeForStrand(refid, '-');
    fw_tree->indexTree();
    rc_tree->indexTree();
  }

  bool
  g2tTree::getGuideExonForTid(int refid, char strand, tid_t tid,
                              uint32_t start, uint32_t end, GuideExon &gexon) {
    IntervalTree *tree = getTreeForStrand(refid, strand);
    if (!tree) return false;
    return tree->findOverlappingForTid(start, end, tid, gexon);
  }

  // Find all guide TIDs that overlap with a read exon
  bool
  g2tTree::getGuideExons(int refid, char strand, GSeg exon,
                        ReadEvaluationConfig config, ExonStatus status,
                        std::vector<GuideExon> &gexons) {
    IntervalTree *tree = getTreeForStrand(refid, strand);
    if (!tree) {
      return false;
    }
    return tree->findOverlapping(exon.start, exon.end, strand, 
      config, status, gexons);
  }
  
};