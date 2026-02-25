
#pragma once

#include "quill/Backend.h"
#include "quill/Frontend.h"
#include "quill/LogMacros.h"
#include "quill/Logger.h"
#include "quill/sinks/ConsoleSink.h"

#define BAM_CMATCH_OVERRIDE 10
#define BAM_CDEL_OVERRIDE 11
#define BAM_CINS_OVERRIDE 12
#define BAM_CLIP_OVERRIDE 13

namespace bramble {

  struct CReadAln;
  struct BamIO;
  struct g2tTree;
  
  struct Cigar {
    uint32_t *cigar = nullptr;
    uint32_t n_cigar = 0;
    uint32_t capacity = 0;

    ~Cigar() {
      delete[] cigar;
    }

    Cigar() = default;

    // copy constructor (deep copy)
    Cigar(const Cigar &other) {
      n_cigar = other.n_cigar;
      capacity = other.capacity;

      if (capacity) {
        cigar = new uint32_t[capacity];
        std::memcpy(cigar, other.cigar, n_cigar * sizeof(uint32_t));
      }
    }

    // copy assignment (deep copy)
    Cigar& operator=(const Cigar &other) {
      if (this == &other) return *this;

      delete[] cigar;

      n_cigar = other.n_cigar;
      capacity = other.capacity;

      cigar = capacity ? new uint32_t[capacity] : nullptr;

      if (n_cigar) {
        std::memcpy(cigar, other.cigar, n_cigar * sizeof(uint32_t));
      }

      return *this;
    }

    // move constructor (cheap, no copy)
    Cigar(Cigar &&other) noexcept {
      cigar = other.cigar;
      n_cigar = other.n_cigar;
      capacity = other.capacity;

      other.cigar = nullptr;
      other.n_cigar = 0;
      other.capacity = 0;
    }

    // move assignment
    Cigar& operator=(Cigar &&other) noexcept {
      if (this == &other) return *this;

      delete[] cigar;

      cigar = other.cigar;
      n_cigar = other.n_cigar;
      capacity = other.capacity;

      other.cigar = nullptr;
      other.n_cigar = 0;
      other.capacity = 0;

      return *this;
    }

    void reserve(uint32_t new_capacity) {
      if (new_capacity <= capacity) return;

      uint32_t *new_buf = new uint32_t[new_capacity];

      if (n_cigar) {
        std::memcpy(new_buf, cigar, n_cigar * sizeof(uint32_t));
      }

      delete[] cigar;
      cigar = new_buf;
      capacity = new_capacity;
    }

    void grow() {
      uint32_t new_capacity = capacity ? capacity * 2 : 4;
      reserve(new_capacity);
    }

    void add_operation(uint32_t len, uint8_t op) {
      if (n_cigar == 0) {
        if (n_cigar == capacity) grow();
        cigar[n_cigar++] = bam_cigar_gen(len, op);
        return;
      }

      uint32_t &prev_cig = cigar[n_cigar - 1];

      uint8_t prev_op = bam_cigar_op(prev_cig);
      uint32_t prev_len = bam_cigar_oplen(prev_cig);

      if (prev_op == op) {
        prev_cig = bam_cigar_gen(prev_len + len, op);
      } else {
        if (n_cigar == capacity) grow();
        cigar[n_cigar++] = bam_cigar_gen(len, op);
      }
    }

  };

  struct GuideExon {
    tid_t tid;
    uint32_t start;
    uint32_t end;
    pos_t    pos;
    uint32_t pos_start;
    uint8_t  exon_id;
    int32_t  left_ins;
    int32_t  right_ins;
    int32_t  left_gap;
    int32_t  right_gap;
    bool     has_prev;
    bool     has_next;
    uint32_t prev_start;
    uint32_t prev_end;
    uint32_t next_start;
    uint32_t next_end;
    uint32_t transcript_len;
    const char* seq;
    uint32_t    seq_len;
  };

  struct AlignInfo {
    pos_t fwpos;
    pos_t rcpos;
    char strand;
    Cigar cigar;
    bool is_reverse;
    bool primary_alignment;
    int clip_score;
    double similarity_score;
    int32_t hit_index;

    AlignInfo()
      : primary_alignment(false),
        clip_score(0) {}
  };

  struct ExonChainMatch {
    tid_t tid;
    bool created;
    AlignInfo align;
    double total_coverage;
    double total_operations;
    int32_t ref_consumed;
    int32_t junc_hits;
    int32_t transcript_len;
    uint8_t prev_op;

    ExonChainMatch()
      : created(false) {}
  };

  enum ExonStatus {
    FIRST_EXON = 0,
    MIDDLE_EXON = 1,
    LAST_EXON = 2,
    ONLY_EXON = 3,
    INS_EXON = 4,
    GAP_EXON = 5,
    LEFTC_EXON = 6,
    RIGHTC_EXON = 7
  };

  struct Segment {
    bool has_gexon;
    bool has_qexon;
    GuideExon gexon;
    GSeg qexon;

    // gap indicated by has_gexon and !has_qexon
    // ins indicated by !has_gexon and has_qexon

    ExonStatus status;
    bool is_small_exon;
    Cigar cigar;
    int score;
  };

  struct TidData {
    bool elim = false;
    bool has_left_clip;
    bool has_right_clip;
    ExonChainMatch match;
    std::vector<Segment> segments;

    TidData()
      : has_left_clip(false),
        has_right_clip(false) {}
  };

  struct ReadOut {
    read_id_t index;
    uint32_t nh;
    uint32_t mapq;
    std::shared_ptr<GSamRecord> brec;

    ReadOut() 
      : index(), 
        nh(), 
        mapq(), 
        brec(nullptr) {}
  };

  struct ReadInfo {
    std::vector<ExonChainMatch> matches;
    bool valid_read;
    bool is_paired;

    std::shared_ptr<ReadOut> read;

    ReadInfo() 
      : matches(), 
        valid_read(false), 
        is_paired(false), 
        read() {}

    ~ReadInfo() {}
  };

  struct BamInfo {
    bool same_transcript;
    bool valid_pair;
    bool is_paired;

    // two ends of a pair (or just one if unpaired)
    std::shared_ptr<ReadOut> read1;
    std::shared_ptr<ReadOut> read2;

    // Read 1 match info
    tid_t r_tid;
    AlignInfo r_align;

    // Read 2 match info
    tid_t m_tid;
    AlignInfo m_align;

    BamInfo() 
      : same_transcript(false), 
        valid_pair(false), 
        is_paired(false) {}

    ~BamInfo() {}
  };

  struct ReadEvaluationConfig {
    uint32_t max_clip;              // max soft clip size
    uint32_t max_ins;               // max insertion to intervals
      // for when there exist no guides to explain a portion of an exon
    uint32_t max_gap;               // max gap in reference to intervals
    bool ignore_small_exons;        // should we ignore small exons?
    uint32_t small_exon_size;       // size of small exon
    uint32_t max_junc_gap;          // max junc mismatch
    float similarity_threshold;     // min similarity score
    //ReadEvaluationResult default_result;
    bool print;
    std::string name;
  };

  struct kswResult {
    uint32_t* cigar_array;      // raw CIGAR ops
    int n_cigar;                // number of CIGAR ops
    int score;
    int max;
  };

  class ReadEvaluator {
  public:
    virtual ~ReadEvaluator();

    // Evaluate a single read and return structured result
    virtual std::vector<ExonChainMatch> 
    evaluate(CReadAln &read, std::shared_ptr<g2tTree> g2t,
            uint8_t *seq, int seq_len) = 0;
  
  protected:
    
    std::vector<char> get_strands_to_check(CReadAln &read);

    void get_clips(CReadAln &read, ReadEvaluationConfig &config, 
                  bool &failure, bool &has_left_clip, bool &has_right_clip,
                  uint32_t &n_left_clip, uint32_t &n_right_clip);

    bool correct_for_gaps(TidData &td, tid_t tid,
                          Segment &seg,
                          ReadEvaluationConfig &config, 
                          std::shared_ptr<g2tTree> g2t,
                          char strand, refid_t refid);

    void get_intervals(unordered_map<tid_t, TidData> &data, CReadAln &read, 
                      uint32_t j, uint32_t exon_count,
                      ReadEvaluationConfig& config, std::shared_ptr<g2tTree> g2t,
                      refid_t refid, char strand, bool has_left_clip,
                      bool has_right_clip, bool& failure);


    void left_clip_rescue(TidData &td, char strand, std::shared_ptr<g2tTree> g2t,
                          refid_t refid, tid_t tid, uint32_t n_left_clip,
                          ReadEvaluationConfig &config, 
                          CReadAln &read, uint8_t *seq, int seq_len);

    void right_clip_rescue(TidData &td, char strand, std::shared_ptr<g2tTree> g2t,
                          refid_t refid, tid_t tid, uint32_t n_right_clip,
                          ReadEvaluationConfig &config, 
                          CReadAln &read, uint8_t *seq, int seq_len);

    void create_match(TidData &td, GuideExon &gexon, tid_t tid, char strand);

    void build_cigar_match(Segment &seg, TidData &td,
                          ExonChainMatch &match, bool first_match,
                          bool last_match, ReadEvaluationConfig &config);

    void build_cigar_ins(Segment &seg, uint32_t k, uint32_t n, 
                        ExonChainMatch &match, ReadEvaluationConfig &config);

    void build_cigar_gap(Segment &seg, ExonChainMatch &match, 
                        ReadEvaluationConfig &config);

    void build_cigar_clip(Segment &seg, ExonChainMatch &match, 
                          ReadEvaluationConfig &config);

    void filter_by_similarity(std::vector<ExonChainMatch> &matches,
                              std::shared_ptr<g2tTree> g2t, 
                              ReadEvaluationConfig config);

  public:
    std::vector<ExonChainMatch>
    evaluate_exon_chains(CReadAln &read, std::shared_ptr<g2tTree> g2t, 
                        ReadEvaluationConfig config,
                        uint8_t *seq, int seq_len);
  };

  class ShortReadEvaluator : public ReadEvaluator {

  public:
    std::vector<ExonChainMatch> 
    evaluate(CReadAln &read, std::shared_ptr<g2tTree> g2t,
            uint8_t *seq, int seq_len) override;

  };

  class LongReadEvaluator : public ReadEvaluator {

  public:
    std::vector<ExonChainMatch> 
    evaluate(CReadAln &read, std::shared_ptr<g2tTree> g2t,
            uint8_t *seq, int seq_len) override;
  };

  void convert_reads(std::vector<CReadAln> &reads,
                    std::shared_ptr<g2tTree> g2t, 
                    std::shared_ptr<ReadEvaluator> evaluator, 
                    quill::Logger *logger,
                    BamIO *io);

}