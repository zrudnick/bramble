
#pragma once

namespace bramble {

  struct CReadAln;
  struct BamIO;
  struct g2tTree;
  struct BundleData;
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
    tid_t tid;
    AlignInfo align;
    pos_t secondary_pos;
    uint32_t transcript_length;
    uint32_t ref_consumed;
    uint32_t total_coverage;
    uint32_t total_read_length;
  };

  struct ReadOut {
    read_id_t index;
    int32_t nh;
    std::string name;
    uint32_t read_size;
    GSamRecord *brec;
    bool is_reverse;
    bool discard_read;

    ReadOut() 
      : index(), 
        nh(), 
        name(),
        read_size(), 
        brec(nullptr), 
        is_reverse(false),  
        discard_read(false) {}

    // brec is deleted by CReadAln
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

  struct ReadEvaluationResult {
    bool valid;
    std::vector<ExonChainMatch> matches;
  };

  struct ReadEvaluationConfig {
    uint32_t max_clip_size;         // max soft clip size
    bool filter_similarity;         // should we filter using similarity score?
    double similarity_threshold;    // similarity threshold
    bool ignore_small_exons;        // should we ignore small exons?
    uint32_t small_exon_size;       // size of small exon
    ReadEvaluationResult default_result;
  };

  class ReadEvaluator {
  public:
    virtual ~ReadEvaluator();

    // Evaluate a single read and return structured result
    virtual ReadEvaluationResult 
    evaluate(BundleData *bundle, uint read_index, g2tTree *g2t) = 0;
  
  protected:
    
    /**
    * Get the match position for this transcript
    *
    * @param interval first interval the read matched to
    * @param tid transcript id
    * @param read_strand strand that the read aligned to
    * @param g2t g2t tree for bundle
    * @param exon_start start of first read exon
    * @param used_backwards_overhang did we match backwards to a previous 
    * interval?
    * (USE_FASTA mode)
    */
    pos_t get_match_pos(IntervalNode *interval, tid_t tid, 
                        char strand, g2tTree *g2t, 
                        uint exon_start, bool used_backwards_overhang);

    std::tuple<uint32_t, uint32_t> 
    get_exon_coordinates(BundleData *bundle, GSeg exon, char strand);

    std::string reverse_complement(const std::string &seq);

    std::string extract_sequence(char *gseq, uint start, 
                                uint length, char strand);

    /**
     * Check if we can match forwards to another interval (USE_FASTA mode)
     *
     * @param interval first interval the read matched to
     */
    bool search_forward(IntervalNode *interval, uint exon_end,
                        char read_strand,
                        std::set<tid_t> &exon_tids, g2tTree *g2t,
                        BundleData *bundle);

    // Function to check backward overhang (left extension)
    bool search_backward(IntervalNode *interval, uint exon_start,
                        char read_strand,
                        std::set<tid_t> &exon_tids, g2tTree *g2t,
                        BundleData *bundle);
    
    std::vector<char> get_strands_to_check(CReadAln* read);

    uint32_t get_boundary_tolerance(uint32_t read_length);

    std::set<tid_t> get_candidate_tids(std::vector<IntervalNode *> intervals,
                                      std::vector<ExonChainMatch> &matches,
                                      std::unordered_map<tid_t, std::shared_ptr<Cigar>> &subcigars,
                                      uint32_t exon_start, uint32_t exon_end,
                                      bool is_first_exon, bool is_last_exon, g2tTree *g2t);

    void add_operation(std::vector<ExonChainMatch> &matches,
                      uint32_t length, uint8_t op, g2tTree *g2t, char strand);

    void hide_small_exon(std::vector<ExonChainMatch> &matches, 
                        uint32_t exon_start, uint32_t exon_end,
                        bool is_first_exon, bool is_last_exon, g2tTree *g2t, char strand);

    bool check_first_exon(uint32_t exon_start, uint32_t interval_start,
                          uint32_t exon_end, uint32_t interval_end,
                          bool is_last_exon, uint32_t max_clip_size);

    void filter_tids(std::vector<IntervalNode *> intervals,
                    std::unordered_map<tid_t, IntervalNode *> prev_intervals,
                    std::set<tid_t> &candidate_tids, 
                    uint32_t exon_start, g2tTree* g2t, 
                    char strand, std::unordered_map<tid_t, uint32_t> &gaps);

    bool check_middle_exon(uint32_t exon_start, uint32_t interval_start,
                          uint32_t exon_end, uint32_t interval_end);

    bool check_last_exon(uint32_t exon_start, uint32_t interval_start,
                        uint32_t exon_end, uint32_t interval_end,
                        uint32_t max_clip_size);

    void build_exon_cigar(std::vector<ExonChainMatch> &matches,
                          std::unordered_map<tid_t, std::shared_ptr<Cigar>> &subcigars,
                          bool is_first_exon, bool is_last_exon,
                          uint32_t exon_start, uint32_t exon_end,
                          uint32_t interval_start, uint32_t interval_end,
                          std::unordered_map<tid_t, uint32_t> &gaps, g2tTree* g2t, char strand);

    void compile_matches(std::vector<ExonChainMatch> &matches,
                        std::set<tid_t> candidate_tids, 
                        std::unordered_map<tid_t, IntervalNode *> &first_intervals, 
                        char strand, g2tTree* g2t, 
                        uint32_t exon_start, uint32_t exon_end,
                        uint32_t cumulative_length, uint32_t coverage_augment,
                        bool used_backwards_overhang);

    bool update_matches(std::vector<ExonChainMatch> &matches,
                        std::set<tid_t> candidate_tids, 
                        uint32_t exon_start, uint32_t exon_end,
                        bool is_first_exon,
                        std::unordered_map<tid_t, uint32_t> &gaps);

    void get_prev_intervals(std::vector<IntervalNode *> intervals,
                            std::set<tid_t> candidate_tids,
                            std::unordered_map<tid_t, IntervalNode *> &prev_intervals);

    void get_first_intervals(std::vector<IntervalNode *> intervals,
                            std::set<tid_t> candidate_tids,
                            std::unordered_map<tid_t, IntervalNode *> &first_intervals);

    void filter_by_similarity(std::vector<ExonChainMatch> &matches,
                              double similarity_threshold);

  public:
    ReadEvaluationResult 
    evaluate_exon_chains(BundleData *bundle, read_id_t id, g2tTree *g2t,
                        ReadEvaluationConfig config);
  };

  class ShortReadEvaluator : public ReadEvaluator {

  public:
    ReadEvaluationResult evaluate(BundleData *bundle, read_id_t id, 
                                  g2tTree *g2t) override;

  };

  class LongReadEvaluator : public ReadEvaluator {

  public:
    ReadEvaluationResult evaluate(BundleData *bundle, read_id_t id, 
                                  g2tTree *g2t) override;
  };

  // -------- function definitions

  void convert_reads(BundleData *bundle, BamIO *io);

}