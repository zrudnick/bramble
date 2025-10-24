
#pragma once

namespace bramble {

  struct BamIO;

  struct IntervalNode {
    uint32_t start, end;
    uint32_t max_end; // maximum end value in subtree
    uint32_t height;
    std::set<tid_t> tids;
    std::unordered_map<tid_t, std::pair<uint32_t, uint32_t>> tid_cum_len;
      // cumulative lengths on left and right sides

    IntervalNode *left;
    IntervalNode *right;

    // Maps from TID to previous/next nodes in genomic order for that TID
    std::unordered_map<tid_t, IntervalNode *>
        tid_prev_nodes; // TID -> previous node with same TID
    std::unordered_map<tid_t, IntervalNode *>
        tid_next_nodes; // TID -> next node with same TID

    IntervalNode(uint32_t s, uint32_t e);
  };

  class IntervalTree {
  public:
    IntervalNode *root;
    IntervalTree();

    ~IntervalTree();

  private:
    // IntervalNode *root;
    std::set<tid_t> all_tids;
    std::unordered_map<tid_t, uint32_t> transcript_lengths;

    int getHeight(IntervalNode *node);

    int getBalance(IntervalNode *node);

    void updateHeight(IntervalNode *node);

    void updateMaxEnd(IntervalNode *node);

    IntervalNode *rotateRight(IntervalNode *y);

    IntervalNode *rotateLeft(IntervalNode *x);

    IntervalNode *getParent(IntervalNode *node);

    IntervalNode *findParent(IntervalNode *current, 
                            IntervalNode *target);

    void destroyTree(IntervalNode *node);

    void inorderTraversal(IntervalNode *node, 
                          std::vector<IntervalNode *> &result);

    IntervalNode *insertNodeBalanced(IntervalNode *node, IntervalNode *newNode);

  public:
    void insert(uint32_t start, uint32_t end, const tid_t &tid);

  private:
    void insertNode(IntervalNode *newNode);

    std::vector<IntervalNode *> findAllOverlapping(uint32_t start, uint32_t end);

    void findOverlappingHelper(IntervalNode *node, uint32_t start, 
                              uint32_t end, std::vector<IntervalNode *> &result);

    std::vector<std::pair<uint32_t, uint32_t>>
    mergeRanges(std::vector<std::pair<uint32_t, uint32_t>> &ranges);

    // Split interval in two
    void splitInterval(IntervalNode *node, uint32_t split_point);

    void getTidNodesHelper(IntervalNode *node, const tid_t &tid,
                          std::vector<IntervalNode *> &result);

    // Get all nodes containing a specific TID
    std::vector<IntervalNode *> getTidNodes(const tid_t &tid);

    // Build chain for a specific TID
    void buildTidChain(const tid_t &tid);

    // Precompute cumulative lengths so that match positions can be calculated
    void precomputeCumulativeLengths(const tid_t &tid, char strand);

  public:
    // For debugging (print_tree)
    std::vector<IntervalNode *> getOrderedIntervals();

    // Find all intervals that overlap with the given range
    std::vector<IntervalNode *> 
    findOverlapping(uint32_t start, uint32_t end, bool allow_gaps);

    // Get next node in chain for a specific TID
    IntervalNode *getNextNodeForTid(IntervalNode *node, const tid_t &tid);

    // Get previous node in chain for a specific TID
    IntervalNode *getPrevNodeForTid(IntervalNode *node, const tid_t &tid);

    // Build all TID chains after tree construction is complete
    void buildAllTidChains();

    // Precompute for all TIDs in the tree
    void precomputeAllCumulativeLengths(char strand);

    // Check for cumulative length
    std::pair<uint32_t, uint32_t>
    findCumulativeLength(IntervalNode *node, const tid_t &tid);

    uint32_t getTranscriptLength(tid_t tid);
  };

  // g2t tree using interval tree
  struct g2tTree {
    // map from transcript names to transcript IDs
    std::unordered_map<tid_t, std::string> name_id_map;
    std::vector<std::string> tid_names;
    const std::string invalid_name{"INVALID TID"};

    IntervalTree *fw_tree; // tree for forward strand guides
    IntervalTree *rc_tree; // tree for reverse strand guides

    g2tTree();

    ~g2tTree();

  private:
    IntervalTree *getTreeForStrand(char strand);

  public:
    // If tid_name is already part of the tree, then return the existing
    // tid for it, otherwise give it the next available tid and return that.
    tid_t insertTidString(const char *&tid_name, BamIO* io);

    // Get the string name given an id
    const std::string &getTidName(tid_t id);

    // Add guide exon with TID and transcript start
    void addGuideExon(uint32_t start, uint32_t end,
                      const tid_t &tid, // const std::string &tid,
                      char strand);

    // Call this after loading all transcript data
    void buildAllTidChains();

    // Call this after loading all transcript data
    void precomputeAllCumulativeLengths();

    // Find all guide TIDs that overlap with a read exon
    std::vector<IntervalNode *> 
    getIntervals(uint32_t readStart, uint32_t readEnd, 
                char strand, bool allow_gaps = false);

    // Get cumulative previous size of exons from transcript
    std::pair<uint32_t, uint32_t> 
    getCumulativeLength(IntervalNode *node, const tid_t &tid, char strand);

    // Get next node in chain for a specific TID
    IntervalNode 
    *getNextNode(IntervalNode *node, const tid_t &tid, char strand);

    // Get previous node in chain for a specific TID
    IntervalNode 
    *getPrevNode(IntervalNode *node, const tid_t &tid, char strand);

    uint32_t getTranscriptLength(const tid_t &tid, char strand);
  };

}