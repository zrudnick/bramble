
#include <numeric>
#include <algorithm>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <vector>

#include <IITree.h>
#include "types.h"
#include "bramble.h"
#include "evaluate.h"
#include "g2t.h"
#include "sw.h"

extern bool SOFT_CLIPS;
extern bool STRICT;
extern bool USE_FASTA;

extern GFastaDb *gfasta;  // FASTA file, -S

namespace bramble {

  IITData::IITData() {}

  GuideExon::GuideExon(uint32_t s, uint32_t e)
    : start(s),
      end(e),
      pos(0),
      pos_start(0),
      exon_id(0),
      left_ins(0),
      right_ins(0),
      left_gap(0),
      right_gap(0),
      seq() {}

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

    // could just grab this later when necessary
    if (USE_FASTA) {
      GFaSeqGet* faseq = gfasta->fetch(ref_name);
      std::string seq = faseq->copyRange(interval.start, 
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

  std::shared_ptr<GuideExon>
  IntervalTree::findOverlappingForTid(uint32_t qstart, uint32_t qend, tid_t tid) {
    
    if (qstart == 0 && qend == 0) return nullptr;
    std::vector<size_t> hits;
    iit->overlap(qstart, qend, hits);

    for (size_t idx : hits) {
      uint32_t s = iit->start(idx);
      uint32_t e = iit->end(idx);

      const IITData& data = iit->data(idx);
      if (data.tid == tid) {
        auto exon = std::make_shared<GuideExon>(s, e);
        exon->transcript_len = data.transcript_len;
        exon->pos_start = data.pos_start;
        exon->exon_id = data.exon_id;
        exon->has_prev = data.has_prev;
        exon->has_next = data.has_next;
        exon->prev_start = data.prev_start;
        exon->prev_end = data.prev_end;
        exon->next_start = data.next_start;
        exon->next_end = data.next_end;
        exon->seq = data.seq; // not sure if this has to be cleared

        return exon;
      }
    }

    return nullptr;
  }

  unordered_map<tid_t, std::shared_ptr<GuideExon>>
  IntervalTree::findOverlapping(uint32_t qstart, uint32_t qend, 
                                char strand, ReadEvaluationConfig config, 
                                ExonStatus status, bool has_left_clip,
                                bool has_right_clip) {
    unordered_map<tid_t, std::shared_ptr<GuideExon>> exons;
    std::vector<size_t> hits;
    
    iit->overlap(qstart, qend, hits);
    if (hits.empty()) {
      if (config.print) {
        printf("=== QUERY FAILED ===\n");
        printf("Query range: [%u, %u], strand: %c\n", qstart, qend, strand);
        printf("Reason: No overlapping intervals found in tree\n\n");
      }
      return exons;
    }

    // Track rejection reasons for debugging
    std::vector<std::string> rejection_reasons;

    for (size_t idx : hits) {
      uint32_t s = iit->start(idx);
      uint32_t e = iit->end(idx);

      if (config.print) printf("tid = %d\n", iit->data(idx).tid);

      // Check soft clip conditions
      if ((!SOFT_CLIPS || STRICT) && (qstart < s || e < qend)) {
        rejection_reasons.push_back(
          "Interval [" + std::to_string(s) + ", " + std::to_string(e) + 
          "] rejected by soft clip policy"
        );
        continue;
      }

      const IITData& data = iit->data(idx);
      auto exon = std::make_shared<GuideExon>(s, e);

      bool rejected = false;
      std::string reject_reason;

      if (strand == '+') {
        if (s <= qstart) {
          // Left junction gap: query starts after interval starts
          exon->pos = (qstart - s) + data.pos_start;
          exon->left_gap = qstart - s;
          if (status == MIDDLE_EXON || status == LAST_EXON) {
            if (exon->left_gap > config.max_junc_gap) {
              reject_reason = "Interval [" + std::to_string(s) + ", " + std::to_string(e) + 
                            "] left junction gap too large: " + std::to_string(exon->left_gap) + 
                            " (max: " + std::to_string(config.max_junc_gap) + ")";
              rejected = true;
            }
          }
        } else {
          exon->pos = data.pos_start;
          // Left clip: query extends before interval
          exon->left_ins = s - qstart;
          if (status == MIDDLE_EXON || status == LAST_EXON) {
            if (exon->left_ins > config.max_ins) {
              reject_reason = "Interval [" + std::to_string(s) + ", " + std::to_string(e) + 
                            "] left insertion too large: " + std::to_string(exon->left_ins) + 
                            " (max: " + std::to_string(config.max_clip) + ")";
              rejected = true;
            }
          } else {
            if (exon->left_ins > config.max_clip) {
              reject_reason = "Interval [" + std::to_string(s) + ", " + std::to_string(e) + 
                            "] left clip too large: " + std::to_string(exon->left_ins) + 
                            " (max: " + std::to_string(config.max_clip) + ")";
              rejected = true;
            }
          }
        }

        if (!rejected) {
          // Check right side
          if (e < qend) {
            // Right clip: query extends past interval
            exon->right_ins = qend - e;
            if (status == FIRST_EXON || status == MIDDLE_EXON) {
              if (exon->right_ins > config.max_ins) {
                reject_reason = "Interval [" + std::to_string(s) + ", " + std::to_string(e) + 
                              "] right insertion too large: " + std::to_string(exon->right_ins) + 
                              " (max: " + std::to_string(config.max_clip) + ")";
                rejected = true;
              }
            } else {
              if (exon->right_ins > config.max_clip) {
                reject_reason = "Interval [" + std::to_string(s) + ", " + std::to_string(e) + 
                              "] right clip too large: " + std::to_string(exon->right_ins) + 
                              " (max: " + std::to_string(config.max_clip) + ")";
                rejected = true;
              }
            }
          } else if (qend < e) {
            // Right junction gap: query ends before interval ends
            exon->right_gap = e - qend;
            if (status == FIRST_EXON || status == MIDDLE_EXON) {
              if (exon->right_gap > config.max_junc_gap) {
                reject_reason = "Interval [" + std::to_string(s) + ", " + std::to_string(e) + 
                              "] right junction gap too large: " + std::to_string(exon->right_gap) + 
                              " (max: " + std::to_string(config.max_junc_gap) + ")";
                rejected = true;
              }
            }
          }
        }

      } else {  // strand == '-'
        if (qend <= e) {
          exon->pos = (e - qend) + data.pos_start;
          exon->right_gap = e - qend;
          // Right junction gap: query ends before interval ends
          if (status == FIRST_EXON || status == MIDDLE_EXON) {
            if (exon->right_gap > config.max_junc_gap) {
              reject_reason = "Interval [" + std::to_string(s) + ", " + std::to_string(e) + 
                            "] right junction gap too large: " + std::to_string(exon->right_gap) + 
                            " (max: " + std::to_string(config.max_junc_gap) + ")";
              rejected = true;
            }
          }
        } else {
          exon->pos = data.pos_start;
          // Right clip: query extends past interval
          exon->right_ins = qend - e;
          if (status == FIRST_EXON || MIDDLE_EXON) {
            if (exon->right_ins > config.max_ins) {
              reject_reason = "Interval [" + std::to_string(s) + ", " + std::to_string(e) + 
                            "] right clip too large: " + std::to_string(exon->right_ins) + 
                            " (max: " + std::to_string(config.max_clip) + ")";
              rejected = true;
            }
          } else {
            if (exon->right_ins > config.max_clip) {
              reject_reason = "Interval [" + std::to_string(s) + ", " + std::to_string(e) + 
                            "] right clip too large: " + std::to_string(exon->right_ins) + 
                            " (max: " + std::to_string(config.max_clip) + ")";
              rejected = true;
            }
          }
        }

        if (!rejected) {
          // Check left side
          if (qstart < s) {
            // Left clip: query extends before interval
            exon->left_ins = s - qstart;
            if (status == MIDDLE_EXON || status == LAST_EXON) {
              if (exon->left_ins > config.max_ins) {
                reject_reason = "Interval [" + std::to_string(s) + ", " + std::to_string(e) + 
                              "] left clip too large: " + std::to_string(exon->left_ins) + 
                              " (max: " + std::to_string(config.max_clip) + ")";
                rejected = true;
              }
            } else {
              if (exon->left_ins > config.max_clip) {
                reject_reason = "Interval [" + std::to_string(s) + ", " + std::to_string(e) + 
                              "] left clip too large: " + std::to_string(exon->left_ins) + 
                              " (max: " + std::to_string(config.max_clip) + ")";
                rejected = true;
              }
            }
            
          } else if (s < qstart) {
            exon->left_gap = qstart - s;
            // Left junction gap: query starts after interval starts
            if (status == MIDDLE_EXON || status == LAST_EXON) {
              if (exon->left_gap > config.max_junc_gap) {
                reject_reason = "Interval [" + std::to_string(s) + ", " + std::to_string(e) + 
                              "] left junction gap too large: " + std::to_string(exon->left_gap) + 
                              " (max: " + std::to_string(config.max_junc_gap) + ")";
                rejected = true;
              }
            }
          }
        }
      }

      if (rejected) {
        rejection_reasons.push_back(reject_reason);
        continue;
      }

      exon->pos_start = data.pos_start;
      exon->exon_id = data.exon_id;

      exon->has_prev = data.has_prev;
      exon->has_next = data.has_next;
      exon->prev_start = data.prev_start;
      exon->prev_end = data.prev_end;
      exon->next_start = data.next_start;
      exon->next_end = data.next_end;

      exon->transcript_len = data.transcript_len;
      exon->seq = data.seq; // not sure if this has to be cleared

      exons[data.tid] = exon;
    }

    // Only print diagnostics if no intervals were found and printing is enabled
    // if (config.print && intervals.empty()) {
    if (config.print) {
      //printf("=== QUERY FAILED ===\n");
      printf("Query range: [%u, %u], strand: %c, status: %d\n", 
             qstart, qend, strand, static_cast<int>(status));
      // printf("Found %zu overlapping intervals in tree, but all were rejected:\n\n", 
      //        hits.size());
      printf("Found %zu overlapping intervals in tree\n\n", 
             hits.size());

      for (size_t idx : hits) {
        uint32_t s = iit->start(idx);
        uint32_t e = iit->end(idx);
        printf("s = %d, e = %d\n", s, e);
      }
      
      for (size_t i = 0; i < rejection_reasons.size(); i++) {
        printf("  %zu. %s\n", i + 1, rejection_reasons[i].c_str());
      }
      printf("\n");
    }

    return exons;
  }

  g2tTree::g2tTree() = default;
    
  g2tTree::~g2tTree() {
    // Clean up all allocated trees
    for (auto& [refid, pair] : trees) {
      delete pair.first;   // forward tree
      delete pair.second;  // reverse tree
    }
    trees.clear();
  }

  // Get or create tree pair for a given refid
  std::pair<IntervalTree*, IntervalTree*>& g2tTree::getTrees(uint8_t refid) {
    auto it = trees.find(refid);
    if (it == trees.end()) {
      // Create new tree pair for this refid
      IntervalTree* fw_tree = new IntervalTree();
      IntervalTree* rc_tree = new IntervalTree();
      trees[refid] = std::make_pair(fw_tree, rc_tree);
      return trees[refid];
    }
    return it->second;
  }

  IntervalTree *g2tTree::getTreeForStrand(uint8_t refid, char strand) {
    auto pair = getTrees(refid);
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
  void g2tTree::indexTrees(uint8_t refid) {
    IntervalTree *fw_tree = getTreeForStrand(refid, '+');
    IntervalTree *rc_tree = getTreeForStrand(refid, '-');
    fw_tree->indexTree();
    rc_tree->indexTree();
  }

  std::shared_ptr<GuideExon> 
  g2tTree::getGuideExonForTid(uint8_t refid, char strand, tid_t tid, 
                              uint32_t start, uint32_t end) {
    IntervalTree *tree = getTreeForStrand(refid, strand);
    return tree->findOverlappingForTid(start, end, tid);
  }

  // Find all guide TIDs that overlap with a read exon
  unordered_map<tid_t, std::shared_ptr<GuideExon>>
  g2tTree::getGuideExons(uint8_t refid, char strand, GSeg exon, 
                        ReadEvaluationConfig config, ExonStatus status,
                        bool has_left_clip, bool has_right_clip) {
    IntervalTree *tree = getTreeForStrand(refid, strand);
    return tree->findOverlapping(exon.start, exon.end, strand, 
      config, status, has_left_clip, has_right_clip);
  }
  
};