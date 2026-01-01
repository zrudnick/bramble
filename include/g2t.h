
#pragma once
#include <IITree.h>
#include "evaluate.h"

namespace bramble {

  struct BamIO;
  struct ReadEvaluationConfig;

  struct IITData {
    tid_t tid;
    uint8_t exon_id;
    uint32_t cum_len;
      // cumulative length on left side, for match position calculation
      
    IITData(tid_t t, uint8_t ei, uint32_t cl);
  };

  struct IntervalNode {
    uint32_t start;
    uint32_t end;
    pos_t pos;
    tid_t tid;
    // std::unordered_map<tid_t, uint32_t> tid_cum_len;
    //   // cumulative length on left side, for match position calculation
    // std::unordered_map<tid_t, uint8_t> tid_exon_id;
    //   // for each tid, this is exon # what?
    uint32_t cum_len;
    uint8_t exon_id;
    int32_t left_clip;
    int32_t right_clip;
      
    IntervalNode(uint32_t s, uint32_t e);
  };

  class IntervalTree {
  public:
    IITree<int, IITData> *iit;

    IntervalTree();

    ~IntervalTree();

  public:
    // Add guide exon with TID and transcript start
    void addInterval(uint32_t start, uint32_t end, const tid_t &tid, 
                    uint8_t ei, uint32_t cl);

    void indexTree();

    // Find intervals that overlap with the given range
    std::vector<std::shared_ptr<IntervalNode>>
    findOverlapping(uint32_t start, uint32_t end, char strand,
                    ReadEvaluationConfig config, ExonStatus status);

  };

  // g2t tree using interval tree
  struct g2tTree {
    // trees for forward and reverse strand
    std::unordered_map<uint8_t, std::pair<IntervalTree*, IntervalTree*>> trees;

    // map from transcript names to transcript IDs
    std::unordered_map<tid_t, std::string> name_id_map;
    std::vector<std::string> tid_names;
    const std::string invalid_name{"INVALID TID"};

    std::unordered_map<refid_t, uint32_t> ref_len_map;

    g2tTree();

    ~g2tTree();

  private:

    std::pair<IntervalTree*, IntervalTree*>& getTrees(uint8_t refid);

    IntervalTree *getTreeForStrand(uint8_t refid, char strand);

  public:
    // If tid_name is already part of the tree, then return the existing
    // tid for it, otherwise give it the next available tid and return that.
    tid_t insertTidString(const char *&tid_name, BamIO* io);

    // Get the string name given an id
    const std::string &getTidName(tid_t id);

    // Add guide exon with TID and transcript start
    void addInterval(refid_t refid, uint32_t start, uint32_t end,
                     const tid_t &tid, uint8_t ei, uint32_t cl, char strand);

    void addRefLen(refid_t refid, uint32_t ref_len);

    // Index tree after guides have been added
    void indexTrees(uint8_t refid);

    // Find all guide TIDs that overlap with a read exon
    std::vector<std::shared_ptr<IntervalNode>> 
    getIntervals(uint8_t refid, char strand, GSeg exon, 
                 ReadEvaluationConfig config, ExonStatus status);

  };

}