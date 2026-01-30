
#pragma once

#include <list>

namespace bramble {

  struct CReadAln;
  struct BamIO;
  struct g2tTree;
  struct GuideExon;

  struct Cigar {
    std::vector<std::pair<uint32_t, uint8_t>> cigar;

    void add_operation(uint32_t length, uint8_t op) {
      if (length == 0) return;
      if (!cigar.empty() && cigar.back().second == op) {
        cigar.back().first += length;
      } else {
        cigar.emplace_back(length, op);
      }
    }

    void prepend_operation(uint32_t length, uint8_t op) {
      if (length == 0) return;
      if (!cigar.empty() && cigar.front().second == op) {
        cigar.front().first += length;
      } else {
        cigar.insert(cigar.begin(), std::make_pair(length, op));
      }
    }

    // for debugging
    std::string to_string() const {
      std::string result;
      for (const auto& [length, op] : cigar) {
        result += std::to_string(length);
        // Convert BAM op code to CIGAR character
        switch(op) {
          case BAM_CMATCH: result += 'M'; break;
          case BAM_CINS: result += 'I'; break;
          case BAM_CDEL: result += 'D'; break;
          case BAM_CREF_SKIP: result += 'N'; break;
          case BAM_CSOFT_CLIP: result += 'S'; break;
          case BAM_CHARD_CLIP: result += 'H'; break;
          case BAM_CPAD: result += 'P'; break;
          case BAM_CEQUAL: result += '='; break;
          case BAM_CDIFF: result += 'X'; break;
          default: result += '?'; break;
        }
      }
      return result.empty() ? "*" : result;
    }

    void reverse() {
      std::reverse(cigar.begin(), cigar.end());
    }

    void clear() {
      cigar.clear();
    }

  };

  struct AlignInfo {
    pos_t fwpos;
    pos_t rcpos;
    char strand;
    std::shared_ptr<Cigar> cigar;
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
    bool is_valid;
    bool has_gexon;
    bool has_qexon;
    std::shared_ptr<GuideExon> gexon;
    GSeg qexon;

    // gap indicated by has_gexon and !has_qexon
    // ins indicated by !has_gexon and has_qexon

    ExonStatus status;
    bool is_small_exon;
    Cigar cigar;
    int score;

    Segment()
      : is_valid(false),
        gexon(nullptr) {}
  };

  struct TidData {
    bool elim = false;
    std::list<Segment> segments;
    ExonChainMatch match;
    bool has_left_ins;
    bool has_right_ins;
    int32_t n_left_ins;
    int32_t n_right_ins;
    bool has_left_clip;
    bool has_right_clip;

    TidData()
      : has_left_ins(false),
        has_right_ins(false),
        n_left_ins(0),
        n_right_ins(0),
        has_left_clip(false),
        has_right_clip(false) {}
  };

  struct ReadOut {
    read_id_t index;
    uint32_t nh;
    uint32_t mapq;
    std::string name;
    uint32_t read_size;
    GSamRecord *brec;
    bool discard_read;

    ReadOut() 
      : index(), 
        nh(), 
        mapq(),
        name(),
        read_size(), 
        brec(nullptr), 
        discard_read(false) {}

    // brec is deleted by CReadAln
  };

  struct ReadInfo {
    std::unordered_map<tid_t, ExonChainMatch> matches;
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

    // Two ends of a pair (or just one if unpaired)
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

  class ReadEvaluator {
  public:
    virtual ~ReadEvaluator();

    // Evaluate a single read and return structured result
    virtual std::unordered_map<tid_t, ExonChainMatch> 
    evaluate(CReadAln * read, read_id_t id, std::shared_ptr<g2tTree> g2t,
            uint8_t *seq, int seq_len) = 0;
  
  protected:
    
    std::string reverse_complement(const std::string &seq);

    std::string extract_sequence(char *gseq, uint start, 
                                uint length, char strand);
    
    std::vector<char> get_strands_to_check(CReadAln * read);

    void get_clips(CReadAln *read, ReadEvaluationConfig &config, bool &failure,
                  bool &has_left_clip, bool &has_right_clip,
                  uint32_t &n_left_clip, uint32_t &n_right_clip);

    void get_intervals(std::unordered_map<tid_t, TidData> &data,
                      CReadAln *read, uint32_t j, uint32_t exon_count,
                      ReadEvaluationConfig &config, std::shared_ptr<g2tTree> g2t,
                      refid_t refid, char strand, bool has_left_clip,
                      bool has_right_clip, Segment &start_ignore,
                      Segment &end_ignore, bool &has_variants, bool &failure);

    void correct_for_gaps(TidData &tid_data, tid_t tid,
                          ReadEvaluationConfig &config, 
                          std::shared_ptr<g2tTree> g2t,
                          char strand, refid_t refid);

    void left_clip_rescue(TidData &tid_data, Segment &start_ignore,
                          char strand, std::shared_ptr<g2tTree> g2t,
                          refid_t refid, tid_t tid, uint32_t n_left_clip,
                          uint32_t n_left_ins, ReadEvaluationConfig &config, 
                          CReadAln *read, uint8_t *seq, int seq_len);

    void right_clip_rescue(TidData &tid_data, Segment &end_ignore,
                          char strand, std::shared_ptr<g2tTree> g2t,
                          refid_t refid, tid_t tid, uint32_t n_right_clip,
                          uint32_t n_right_ins, ReadEvaluationConfig &config, 
                          CReadAln *read, uint8_t *seq, int seq_len);

    void create_match(std::unordered_map<tid_t, TidData> &data, 
                                  std::shared_ptr<GuideExon> gexon, tid_t tid,
                                  char strand, char read_strand, bool has_gexon,
                                  uint32_t ins_exon_len, ReadEvaluationConfig config);

    void build_cigar_match(Segment &segment, TidData &tid_data,
                          ExonChainMatch &match, bool first_match,
                          bool last_match, ReadEvaluationConfig &config);

    void build_cigar_ins(Segment &segment, TidData &tid_data, uint32_t k,
                        uint32_t n, ExonChainMatch &match, ReadEvaluationConfig &config);

    void build_cigar_gap(Segment &segment, TidData &tid_data,
                        ExonChainMatch &match, ReadEvaluationConfig &config);

    void build_cigar_clip(Segment &segment, TidData &tid_data,
                          ExonChainMatch &match, ReadEvaluationConfig &config);

    void filter_by_similarity(std::unordered_map<tid_t, ExonChainMatch> &matches,
                              std::shared_ptr<g2tTree> g2t, ReadEvaluationConfig config);

  public:
    std::unordered_map<tid_t, ExonChainMatch> 
    evaluate_exon_chains(CReadAln *read, read_id_t id, 
                        std::shared_ptr<g2tTree> g2t, ReadEvaluationConfig config,
                        uint8_t *seq, int seq_len);
  };

  class ShortReadEvaluator : public ReadEvaluator {

  public:
    std::unordered_map<tid_t, ExonChainMatch> 
    evaluate(CReadAln *read, read_id_t id, std::shared_ptr<g2tTree> g2t,
            uint8_t *seq, int seq_len) override;

  };

  class LongReadEvaluator : public ReadEvaluator {

  public:
    std::unordered_map<tid_t, ExonChainMatch> 
    evaluate(CReadAln *read, read_id_t id, std::shared_ptr<g2tTree> g2t,
            uint8_t *seq, int seq_len) override;
  };

  // -------- function definitions

  void convert_reads(std::vector<CReadAln *> &reads,
                    std::shared_ptr<g2tTree> g2t, 
                    std::shared_ptr<ReadEvaluator> evaluator, 
                    BamIO *io);

}