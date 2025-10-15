
#pragma once

namespace bramble {

  struct CReadAln;
  struct BamIO;
  struct g2tTree;
  struct BundleData;
  struct IntervalNode;

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

    void clear() {
      cigar.clear();
    }

  };

  struct ExonChainMatch {
    tid_t tid;
    pos_t pos;
    std::vector<IntervalNode*> matched_intervals;
    uint32_t total_coverage;
    uint32_t gaps_count;
    uint32_t total_gap_size;
    double similarity_score;
    uint32_t soft_clip_front;
    uint32_t soft_clip_back;
  };

  struct ReadOut {
    read_id_t index;
    int32_t nh;
    std::string name;
    uint32_t read_size;
    GSamRecord *brec;
    Cigar cigar;
    bool is_reverse;
    bool discard_read;
    char strand;

    ReadOut() 
      : index(), 
        nh(), 
        name(),
        read_size(), 
        brec(nullptr),
        cigar(), 
        is_reverse(false),  
        discard_read(false),
        strand('.') {}

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
    uint32_t r_pos;

    // Read 2 match info
    tid_t m_tid;
    uint32_t m_pos;

    BamInfo() 
      : same_transcript(false), 
        valid_pair(false), 
        is_paired(false) {}

    ~BamInfo() {}
  };

  struct ReadEvaluationResult {
    bool valid = false;
    std::vector<ExonChainMatch> matches;
    char strand = '.';
    Cigar cigar;
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
                        char read_strand, g2tTree *g2t, 
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

    uint32_t boundary_tolerance(uint32_t read_length);

    std::set<tid_t> get_candidate_tids(std::vector<IntervalNode *> intervals);

    void hide_small_exon(Cigar& cigar, uint32_t exon_start, uint32_t exon_end,
                         bool is_first_exon, bool is_last_exon);

    double 
    exon_chain_similarity(GVec<GSeg>& read_exons,
                          const std::vector<std::vector<IntervalNode*>>& 
                          all_intervals,
                          tid_t tid, char strand, BundleData* bundle);

    bool check_first_exon(uint32_t exon_start, uint32_t interval_start,
                          uint32_t exon_end, uint32_t interval_end,
                          bool is_last_exon, uint32_t max_clip_size);

    void filter_tids(IntervalNode* first_interval, 
                    IntervalNode* prev_last_interval,
                    std::set<tid_t> &candidate_tids, 
                    uint32_t exon_start, g2tTree* g2t, char strand);

    bool check_middle_exon(uint32_t exon_start, uint32_t interval_start,
                          uint32_t exon_end, uint32_t interval_end,
                          IntervalNode* first_interval, 
                          IntervalNode* prev_last_interval,
                          std::set<tid_t> &candidate_tids, g2tTree* g2t, 
                          char strand);

    bool check_last_exon(uint32_t exon_start, uint32_t interval_start,
                        uint32_t exon_end, uint32_t interval_end,
                        IntervalNode* first_interval,
                        IntervalNode* prev_last_interval, 
                        std::set<tid_t> &candidate_tids,
                        g2tTree* g2t, char strand, uint32_t max_clip_size);

    void build_exon_cigar(Cigar& cigar, bool is_first_exon, bool is_last_exon,
                          uint32_t exon_start, uint32_t exon_end,
                          uint32_t interval_start, uint32_t interval_end);

    bool compile_matches(std::vector<ExonChainMatch> &matches,
                        std::set<tid_t> candidate_tids, 
                        IntervalNode *first_interval, 
                        char strand, g2tTree* g2t, uint32_t exon_start, 
                        bool is_first_exon, bool used_backwards_overhang);

    void filter_by_similarity(std::vector<ExonChainMatch> &matches,
                              std::vector<std::vector<IntervalNode*>> 
                              all_matched_intervals,
                              GVec<GSeg> read_exons, char strand, 
                              BundleData *bundle, 
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