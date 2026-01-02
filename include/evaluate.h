
#pragma once

namespace bramble {

  struct CReadAln;
  struct BamIO;
  struct g2tTree;
  struct IntervalNode;

  // i'll make this a uint32_t* at some point
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

    void clear() {
      cigar.clear();
    }

  };

  struct AlignInfo {
    pos_t pos;
    char strand;
    std::shared_ptr<Cigar> cigar;
    bool primary_alignment;
    double similarity_score;
    int32_t hit_index;

    AlignInfo()
      : primary_alignment(false) {}
  };

  struct ExonChainMatch {
    AlignInfo align;
    pos_t secondary_pos;
    uint32_t transcript_length;
    uint32_t ref_consumed;
    uint32_t total_coverage;
    uint32_t total_read_length;
  };

  enum ExonStatus {
    FIRST_EXON = 0,
    MIDDLE_EXON = 1,
    LAST_EXON = 2,
    ONLY_EXON = 3
  };

  struct ReadOut {
    read_id_t index;
    uint32_t nh;
    uint32_t mapq;
    std::string name;
    uint32_t read_size;
    GSamRecord *brec;
    bool is_reverse;
    bool discard_read;

    ReadOut() 
      : index(), 
        nh(), 
        mapq(),
        name(),
        read_size(), 
        brec(nullptr), 
        is_reverse(false),  
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

  // struct ReadEvaluationResult {
  //   bool valid;
  //   std::unordered_map<tid_t, ExonChainMatch> matches;
  // };

  struct ReadEvaluationConfig {
    uint32_t boundary_tolerance;    // boundary tolerance
    uint32_t max_clip_size;         // max soft clip size
    uint32_t max_ins;               // max insertion to intervals
      // for when there exist no guides to explain a portion of an exon
    uint32_t max_gap;               // max gap in reference to intervals
    bool filter_similarity;         // should we filter using similarity score?
    double similarity_threshold;    // similarity threshold
    bool ignore_small_exons;        // should we ignore small exons?
    uint32_t small_exon_size;       // size of small exon
    uint32_t max_junc_gap;          // max junc mismatch
    //ReadEvaluationResult default_result;
    bool print;
  };

  class ReadEvaluator {
  public:
    virtual ~ReadEvaluator();

    // Evaluate a single read and return structured result
    virtual std::unordered_map<tid_t, ExonChainMatch> 
    evaluate(CReadAln * read, read_id_t id, g2tTree *g2t) = 0;
  
  protected:
    
    std::string reverse_complement(const std::string &seq);

    std::string extract_sequence(char *gseq, uint start, 
                                uint length, char strand);
    
    std::vector<char> get_strands_to_check(CReadAln * read);

    void add_operation_for_all(std::unordered_map<tid_t, ExonChainMatch> &matches,
                              uint32_t length, uint8_t op, g2tTree *g2t, 
                              char strand);

    void hide_small_exon(std::unordered_map<tid_t, ExonChainMatch> &matches,
                        uint32_t exon_start, uint32_t exon_end,
                        ExonStatus status, g2tTree *g2t, char strand);

    void ensure_continuity(std::vector<std::shared_ptr<IntervalNode>> &intervals,
                          std::unordered_map<tid_t, uint8_t> &exon_id_map, 
                          ExonStatus status, char strand);

    // bool check_first_exon(std::vector<std::shared_ptr<IntervalNode>> intervals,
    //                       GSeg exon, bool is_last_exon, uint32_t max_clip_size,
    //                       uint32_t tolerance, uint32_t max_junc_gap,
    //                       char strand);

    // bool check_middle_exon(uint32_t exon_start, uint32_t interval_start,
    //                       uint32_t exon_end, uint32_t interval_end,
    //                       uint32_t tolerance, uint32_t max_junc_gap,
    //                       char strand);

    // bool check_last_exon(uint32_t exon_start, uint32_t interval_start,
    //                     uint32_t exon_end, uint32_t interval_end,
    //                     uint32_t max_clip_size, uint32_t tolerance, 
    //                     uint32_t max_junc_gap, char strand);

    void build_exon_cigar(std::unordered_map<tid_t, ExonChainMatch> &matches,
                          ExonStatus status, uint32_t exon_start, uint32_t exon_end,
                          std::vector<std::shared_ptr<IntervalNode>> &intervals,
                          std::unordered_map<tid_t, uint32_t> &gaps, g2tTree* g2t, char strand);

    void compile_matches(std::unordered_map<tid_t, ExonChainMatch> &matches,
                        std::vector<std::shared_ptr<IntervalNode>> &intervals, 
                        char strand, uint32_t exon_start, uint32_t exon_end,
                        uint32_t cumulative_length, uint32_t coverage_augment);

    bool update_matches(std::unordered_map<tid_t, ExonChainMatch> &matches,
                        std::vector<std::shared_ptr<IntervalNode>> &intervals, 
                        char strand, uint32_t exon_start, uint32_t exon_end,
                        ExonStatus status);

    // void get_prev_intervals(std::vector<IntervalNode> intervals,
    //                         std::set<tid_t> candidate_tids,
    //                         std::unordered_map<tid_t, IntervalNode> &prev_intervals,
    //                         char strand);

    // void get_first_intervals(std::vector<IntervalNode> intervals,
    //                         std::set<tid_t> candidate_tids,
    //                         char strand,
    //                         std::unordered_map<tid_t, IntervalNode> &first_intervals);

    void filter_by_similarity(std::unordered_map<tid_t, ExonChainMatch> &matches,
                              double similarity_threshold);

  public:
    std::unordered_map<tid_t, ExonChainMatch> 
    evaluate_exon_chains(CReadAln *read, read_id_t id, 
                        g2tTree *g2t, ReadEvaluationConfig config);
  };

  class ShortReadEvaluator : public ReadEvaluator {

  public:
    std::unordered_map<tid_t, ExonChainMatch> 
    evaluate(CReadAln *read, read_id_t id, g2tTree *g2t) override;

  };

  class LongReadEvaluator : public ReadEvaluator {

  public:
    std::unordered_map<tid_t, ExonChainMatch> 
    evaluate(CReadAln *read, read_id_t id, g2tTree *g2t) override;
  };

  // -------- function definitions

  void convert_reads(std::vector<CReadAln *> &reads,
                    g2tTree *g2t, ReadEvaluator *evaluator, 
                    BamIO *io);

}