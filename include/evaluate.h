
#pragma once

namespace bramble {

  struct CReadAln;
  struct BamIO;
  struct g2tTree;
  struct BundleData;
  struct IntervalNode;

  struct match_hash {
    inline std::size_t operator()(const std::pair<read_id_t, pos_t> & v) const {
        return v.first*31 + v.second;
    }
  };

  struct ReadOut {
    read_id_t index;
    int32_t nh;
    uint32_t read_size;
    GSamRecord *brec;
    bool has_introns;
    bool is_reverse;
    uint32_t soft_clip_front;
    uint32_t soft_clip_back;
    bool discard_read;
    char strand;

    ReadOut() 
      : index(), 
        nh(), 
        read_size(), 
        brec(nullptr), 
        has_introns(false),
        is_reverse(false), 
        soft_clip_front(0), 
        soft_clip_back(0), 
        discard_read(false),
        strand('.') {}
  };

  struct ReadInfo {
    std::unordered_set<std::pair<tid_t, pos_t>, match_hash> matches;
    bool valid_read;
    bool is_paired;

    ReadOut *read;

    ReadInfo() 
      : matches(), 
        valid_read(false), 
        is_paired(false), 
        read() {}

    ~ReadInfo() {
      // brec is deleted by CReadAln
    }
  };

  struct BamInfo {
    bool same_transcript;
    bool valid_pair;
    bool is_paired;

    // Two ends of a pair (or just one if unpaired)
    ReadOut *read1;
    ReadOut *read2;

    tid_t r_tid;
    uint32_t r_pos;
    tid_t m_tid;
    uint32_t m_pos;

    BamInfo() : same_transcript(false), valid_pair(false), is_paired(false) {}
  };

  struct GroupKey {
    uint start;
    std::string cigar_str;

    // for std::map
    bool operator<(const GroupKey &other) const {
      if (start != other.start)
        return start < other.start;
      return cigar_str < other.cigar_str;
    }

    bool operator==(const GroupKey &other) const {
      return start == other.start && cigar_str == other.cigar_str;
    }
  };

  struct AlnGroups {
    std::map<GroupKey, uint> key_to_group;
    std::vector<std::vector<uint32_t>> groups;

    AlnGroups() 
      : key_to_group(), 
        groups() {}

    std::string getCigar(bam1_t *b) {
      uint32_t n_cigar = b->core.n_cigar;
      uint32_t *cigar = bam_get_cigar(b);

      std::string cigar_str;
      cigar_str.reserve(n_cigar * 8);

      for (uint32_t i = 0; i < n_cigar; i++) {
        uint32_t op_len = bam_cigar_oplen(cigar[i]);
        char op_char = bam_cigar_opchr(cigar[i]);
        cigar_str += std::to_string(op_len) + op_char;
      }
      return cigar_str;
    }

    void Add(CReadAln *read, uint n) {
      if (!read || !read->brec)
        return;
      bam1_t *b = read->brec->get_b();

      std::string cigar = getCigar(b);
      GroupKey key{read->start, cigar};

      auto it = key_to_group.find(key);
      uint32_t group_num;

      // Key exists
      if (it != key_to_group.end()) {
        group_num = it->second;
        groups[group_num].push_back(n);

        // New key, create new group
      } else {
        group_num = groups.size();
        key_to_group[key] = group_num;
        groups.push_back(std::vector<uint>{n});
      }
    }

    void Clear() {
      key_to_group.clear();
      groups.clear();
    }

    ~AlnGroups() { Clear(); }
  };

  struct ReadEvaluationResult {
    bool valid = false;
    std::unordered_set<std::pair<tid_t, pos_t>, match_hash> matches;
    char strand = '.';
    uint32_t soft_clip_front = 0;
    uint32_t soft_clip_back = 0;

    // debug / metrics
    uint32_t longest_overhang = 0;
    std::string reject_reason; // optional
  };

  class ReadEvaluator {
  public:
    virtual ~ReadEvaluator();

    // Evaluate a single read (or read group) and return structured result.
    virtual ReadEvaluationResult 
    evaluate(BundleData *bundle, uint read_index, g2tTree *g2t,
            const std::vector<uint32_t> &group) = 0;
  
  protected:
    
    /**
    * Get the match position for this transcript
    *
    * @param interval first interval the read matched to
    * @param tid transcript id
    * @param read_strand strand that the read aligned to
    * @param g2t g2t tree for bundle
    * @param exon_start start of first read exon
    * @param used_backwards_overhang did we match backwards to a previous interval?
    * (USE_FASTA mode)
    */
    uint get_match_pos(IntervalNode *interval, tid_t tid, char read_strand,
                      g2tTree *g2t, uint exon_start, bool used_backwards_overhang);

    std::string reverse_complement(const std::string &seq);

    std::string extract_sequence(char *gseq, uint start, uint length, char strand);

    /**
     * Check if we can match forwards to another interval (USE_FASTA mode)
     *
     * @param interval first interval the read matched to
     */
    bool check_forward_overhang(IntervalNode *interval, uint exon_end,
                                char read_strand,
                                std::set<tid_t> &exon_tids, g2tTree *g2t,
                                BundleData *bundle);

    // Function to check backward overhang (left extension)
    bool check_backward_overhang(IntervalNode *interval, uint exon_start,
                                char read_strand,
                                std::set<tid_t> &exon_tids, g2tTree *g2t,
                                BundleData *bundle);
  };

  class ShortReadEvaluator : public ReadEvaluator {
  private:
    std::set<tid_t> 
    collapse_intervals(std::vector<IntervalNode *> sorted_intervals,
                      uint exon_start, bool is_first_exon, 
                      char read_strand, g2tTree *g2t,
                      BundleData *bundle, bool &used_backwards_overhang,
                      uint &soft_clip, IntervalNode *prev_last_interval);

  public:
    ReadEvaluationResult evaluate(BundleData *bundle,
                                  uint read_index,
                                  g2tTree *g2t,
                                  const std::vector<uint32_t> &group) override;

  };

  class LongReadEvaluator : public ReadEvaluator {
  public:
    ReadEvaluationResult evaluate(BundleData *bundle,
                                  uint read_index,
                                  g2tTree *g2t,
                                  const std::vector<uint32_t> &group) override;

  private:
    // helper methods: collapse_intervals_with_score, splice_rescue, compute_exon_chain_score, ...
  };

  // -------- function definitions

  void print_tree(g2tTree* g2t);

  std::unique_ptr<g2tTree> make_g2t_tree(BundleData *bundle, BamIO *io);

  void process_read_out(BundleData *&bundle, uint read_index, g2tTree *g2t,
                        std::unordered_map<read_id_t, ReadInfo *> &read_info,
                        std::vector<read_id_t> group,
                        std::unique_ptr<ReadEvaluator> &evaluator);

  void add_mate_info(const std::unordered_set<tid_t> &final_transcripts,
                    const std::unordered_set<tid_t> &read_transcripts,
                    const std::unordered_set<tid_t> &mate_transcripts,
                    const std::unordered_map<tid_t, pos_t> &read_positions,
                    const std::unordered_map<tid_t, pos_t> &mate_positions,
                    std::unordered_map<read_id_t, ReadInfo *> &read_info, 
                    std::unordered_map<bam_id_t, BamInfo *> &bam_info, 
                    uint32_t read_index, uint32_t mate_index, uint32_t mate_case);

  inline uint64_t make_pair_key(read_id_t a, read_id_t b);

  void update_read_matches(ReadInfo *read_info,
                          const std::unordered_set<tid_t> &final_transcripts);

  void process_mate_pairs(BundleData *bundle, 
                          std::unordered_map<read_id_t, ReadInfo *> &read_info,
                          std::unordered_map<bam_id_t, BamInfo *> &bam_info);

  void free_read_data(std::unordered_map<read_id_t, ReadInfo *> &read_info,
                      std::unordered_map<bam_id_t, BamInfo *> &bam_info);

  void convert_reads(BundleData *bundle, BamIO *io);

}