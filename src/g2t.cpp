
#include <numeric>
#include <algorithm>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <IITree.h>
#include "types.h"
#include "bundles.h"
#include "bramble.h"
#include "evaluate.h"
#include "g2t.h"

extern bool SOFT_CLIPS;
extern bool STRICT;

namespace bramble {

  IITData::IITData(tid_t t, uint8_t ei, uint32_t cl)
    : tid(t), 
      exon_id(ei),
      cum_len(cl) {}

  IntervalNode::IntervalNode(uint32_t s = 0, uint32_t e = 0)
    : start(s),
      end(e),
      pos(0),
      tid(0),
      cum_len(0),
      exon_id(0),
      left_clip(0),
      right_clip(0) {}

  IntervalTree::IntervalTree() {
    iit = new IITree<int, IITData>();
  }

  IntervalTree::~IntervalTree() { 
    delete(iit); 
  }

  // Add guide exon with TID and transcript start
  void IntervalTree::addInterval(uint32_t start, uint32_t end,
                                const tid_t &tid, uint8_t ei, 
                                uint32_t cl) {
    IITData node = IITData(tid, ei, cl);
    iit->add(start, end, node);
  }

  void IntervalTree::indexTree() {
    iit->index();
  }

  struct iitData {
    uint32_t start;
    uint32_t end;
    IITData data;
  };

  std::vector<std::shared_ptr<IntervalNode>>
  IntervalTree::findOverlapping(uint32_t qstart, uint32_t qend, 
                                char strand, ReadEvaluationConfig config, 
                                ExonStatus status) {
    
    std::vector<std::shared_ptr<IntervalNode>> intervals;
    std::vector<size_t> hits;
    iit->overlap(qstart, qend, hits);
    if (hits.empty()) {
      // printf("hits are empty\n");
      return intervals;
    }

    for (size_t idx : hits) {
      uint32_t s = iit->start(idx);
      uint32_t e = iit->end(idx);

      if (qstart < s && e < qend) continue;
      if ((!SOFT_CLIPS || STRICT) && (qstart < s || e < qend)) continue;
      const IITData& data = iit->data(idx);

      auto interval = std::make_shared<IntervalNode>(s, e);
      interval->tid = data.tid;

      // if (interval->tid == 14172) {
      //   printf("New exon\n");
      //   printf("s = %d, e = %d\n", s, e);
      //   printf("qs = %d, qe = %d\n", qstart, qend);
      // }

      if (strand == '+') {
        if (s < qstart) {
          interval->pos = (qstart - s) + data.cum_len;
          if (e < qend) {
            interval->right_clip = qend - e;
            if (interval->right_clip > config.max_clip_size) {
              continue;
            }
          }
        } else {
          interval->pos = data.cum_len;
          interval->left_clip = s - qstart;
          if (interval->left_clip > config.max_clip_size) {
            continue;
          }
        }

      } else {
        if (qend < e) {
          interval->pos = (e - qend) + data.cum_len;
          //interval->pos = ()
          // if (interval->tid == 14172) {
          //   printf("qend < e\n");
          //   printf("pos = %d, cum_len = %d\n", interval->pos, data.cum_len);
          // }
          if (qstart < s) {
            interval->left_clip = s - qstart;
            if (interval->left_clip > config.max_clip_size) {
              continue;
            }
          }
        } else {
          interval->pos = data.cum_len;
          // if (interval->tid == 14172) {
          //   printf("e <= qend\n");
          //   printf("pos = %d, cum_len = %d\n", interval->pos, data.cum_len);
          // }
          interval->right_clip = qend - e;
          if (interval->right_clip > config.max_clip_size) {
            continue;
          }
        }
      }

      // Ensure correct splice junctions
      if ((status == FIRST_EXON || status == MIDDLE_EXON) 
          && (qend < e)) { // right side
        if (e - qend > config.max_junc_gap) {
          continue;
        }
      }
      if ((status == MIDDLE_EXON || status == LAST_EXON) 
          && (s < qstart)) { // left side
        if (qstart - e > config.max_junc_gap) {
          continue;
        }
      }

      interval->cum_len = data.cum_len;
      interval->exon_id = data.exon_id;

      intervals.emplace_back(interval);
    }

    //printf("intervals are compiled!\n");
    // for (auto &interval : intervals) {
    //   if (interval->tid == 14172) {
    //     printf("tid = %d\n", interval->tid);
    //   }
    // }

    return intervals;

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
  void g2tTree::addInterval(refid_t refid, uint32_t start, uint32_t end,
                            const tid_t &tid, uint8_t ei, uint32_t cl, char strand) {
    IntervalTree *tree = getTreeForStrand(refid, strand);
    IITData node = IITData(tid, ei, cl);
    tree->addInterval(start, end, tid, ei, cl);
  }

  void g2tTree::addRefLen(refid_t refid, uint32_t ref_len) {
    ref_len_map[refid] = ref_len;
  }

  // Index tree after guides have been added
  void g2tTree::indexTrees(uint8_t refid) {
    IntervalTree *fw_tree = getTreeForStrand(refid, '+');
    IntervalTree *rc_tree = getTreeForStrand(refid, '-');
    fw_tree->indexTree();
    rc_tree->indexTree();
  }

  // Find all guide TIDs that overlap with a read exon
  std::vector<std::shared_ptr<IntervalNode>>
  g2tTree::getIntervals(uint8_t refid, char strand, GSeg exon, 
                        ReadEvaluationConfig config, ExonStatus status) {
    IntervalTree *tree = getTreeForStrand(refid, strand);
    return tree->findOverlapping(exon.start, exon.end, strand, 
      config, status);
  }
  
};