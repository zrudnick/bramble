
#pragma once
#include <IITree.h>
#include "evaluate.h"

namespace bramble {

  struct BamIO;
  struct ReadEvaluationConfig;
  
  struct IntervalData {
    uint32_t start;
    uint32_t end;
    uint8_t idx;
    uint32_t pos_start;

    bool has_prev;
    bool has_next;
    uint32_t prev_start;
    uint32_t prev_end;
    uint32_t next_start;
    uint32_t next_end;

    uint32_t transcript_len;
  };

  struct IITData {
    tid_t tid;
    uint8_t exon_id;
    uint32_t pos_start;
      // cumulative length on left side, for match position calculation
    bool has_prev;
    bool has_next;
    uint32_t prev_start;
    uint32_t prev_end;
    uint32_t next_start;
    uint32_t next_end;
    uint32_t transcript_len;
    std::string seq;
      
    IITData();
  };

  struct GuideExon {
    uint32_t start;
    uint32_t end;
    pos_t pos;  // calculated using alignment info
    uint32_t pos_start;
    uint8_t exon_id;
    int32_t left_ins;
    int32_t right_ins;
    int32_t left_gap;
    int32_t right_gap;

    bool has_prev;
    bool has_next;
    uint32_t prev_start;
    uint32_t prev_end;
    uint32_t next_start;
    uint32_t next_end;

    uint32_t transcript_len;
    std::string seq;
      
    GuideExon(uint32_t s, uint32_t e);
  };

  class IntervalTree {
  public:
    IITree<int, IITData> *iit;

    IntervalTree();

    ~IntervalTree();

  public:
    // Add guide exon with TID and transcript start
    void addInterval(const tid_t &tid, IntervalData interval, 
                    const char* ref_name);

    void indexTree();

    std::shared_ptr<GuideExon>
    findOverlappingForTid(uint32_t start, uint32_t end, tid_t tid);

    // Find intervals that overlap with the given range
    std::unordered_map<tid_t, std::shared_ptr<GuideExon>>
    findOverlapping(uint32_t start, uint32_t end, char strand,
                    ReadEvaluationConfig config, ExonStatus status,
                    bool has_left_clip, bool has_right_clip);

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
    void addInterval(refid_t refid, const tid_t &tid, 
                    IntervalData interval, char strand,
                    const char* ref_name);

    // Index tree after guides have been added
    void indexTrees(uint8_t refid);

    std::shared_ptr<GuideExon> 
    getGuideExonForTid(uint8_t refid, char strand, tid_t tid, 
                      uint32_t start, uint32_t end);

    // Find all guide exons that overlap with a read exon
    std::unordered_map<tid_t, std::shared_ptr<GuideExon>> 
    getGuideExons(uint8_t refid, char strand, GSeg exon, 
                 ReadEvaluationConfig config, ExonStatus status,
                 bool has_left_clip, bool has_right_clip);

  };

}