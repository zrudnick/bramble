#ifndef MAIN_H
#define MAIN_H

#include <algorithm>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <unordered_map>
#include <vector>

#include "GArgs.h"
#include "GBitVec.h"
#include "GHashMap.hh"
#include "GSam.h"
#include "GStr.h"
#include "gff.h"
#include "time.h"

using namespace std;

#define MAX_NODE 1000000
#define SMALL_EXON                                                             \
  35 // exons smaller than this have a tendency to be missed by long read data
#define RUNOFF_DIST 200  // Reads at what distance should be considered part of separate bundles?

#define LONG_INTRON_ANCHOR 25
#define MISMATCH_FRAC 0.02

using tid_t = uint32_t;
using read_id_t = uint32_t;
using bam_id_t = uint32_t;

struct BamIO {
protected:
  std::unique_ptr<GSamReader> reader{nullptr};
  std::unique_ptr<GSamWriter> writer{nullptr};
  GSamRecord *current_record = nullptr;
  GStr input_bam;
  GStr output_bam;
  GStr header_sam;

public:
  BamIO(const GStr &in_bam, const GStr &out_bam, const GStr &sam_header)
      : input_bam(in_bam), output_bam(out_bam), header_sam(sam_header) {}

  void start() {
    reader = std::make_unique<GSamReader>(input_bam.chars());
    writer = std::make_unique<GSamWriter>(output_bam.chars(), header_sam.chars());

    GSamRecord *next_rec = reader->next();
    if (next_rec) {
      current_record = next_rec;
    }
  }

  GSamRecord *next() {
    GSamRecord *result = current_record;

    current_record = reader->next(); // null at EOF

    return result;
  }

  void write(GSamRecord *rec) {
    if (writer && rec != nullptr) {
      writer->write(rec);
    }
  }

  int32_t get_tid(const char *transcript_id) {
    return writer->get_tid(transcript_id);
  }

  void stop() {
    /*
    if (reader != nullptr) {
      delete reader;
      reader = nullptr;
    }
    if (writer != nullptr) {
      delete writer;
      writer = nullptr;
    }
    */
    if (current_record != nullptr) {
      delete current_record;
      current_record = nullptr;
    }
  }
};

struct IntervalNode {
  uint32_t start, end;
  uint32_t max_end; // Maximum end value in subtree (interval tree augmentation)
  uint32_t height;
  std::set<tid_t> tids;
  std::unordered_map<tid_t, uint32_t> tid_cum_len;

  IntervalNode *left;
  IntervalNode *right;

  // Maps from TID to previous/next nodes in genomic order for that TID
  std::unordered_map<tid_t, IntervalNode *>
      tid_prev_nodes; // TID -> previous node with same TID
  std::unordered_map<tid_t, IntervalNode *>
      tid_next_nodes; // TID -> next node with same TID

  IntervalNode(uint32_t s, uint32_t e)
      : start(s), end(e), max_end(e), height(1), left(nullptr), right(nullptr) {
  }
};

class IntervalTree {
private:
  IntervalNode *root;
  std::set<tid_t> all_tids;
  //std::unordered_map<tid_t, std::vector<IntervalNode*>> tid_chains;

  int getHeight(IntervalNode *node) { return node ? node->height : 0; }

  int getBalance(IntervalNode *node) {
    return node ? getHeight(node->left) - getHeight(node->right) : 0;
  }

  void updateHeight(IntervalNode *node) {
    if (node) {
      uint32_t left_height = getHeight(node->left);
      uint32_t right_height = getHeight(node->right);
      node->height = 1 + std::max(left_height, right_height);
    }
  }

  void updateMaxEnd(IntervalNode *node) {
    if (!node)
      return;

    node->max_end = node->end;
    if (node->left) {
      node->max_end = std::max(node->max_end, node->left->max_end);
    }
    if (node->right) {
      node->max_end = std::max(node->max_end, node->right->max_end);
    }
  }

  IntervalNode *rotateRight(IntervalNode *y) {
    IntervalNode *x = y->left;
    IntervalNode *T2 = x->right;

    // Rotation
    x->right = y;
    y->left = T2;

    updateHeight(y);
    updateHeight(x);
    updateMaxEnd(y);
    updateMaxEnd(x);

    return x;
  }

  IntervalNode *rotateLeft(IntervalNode *x) {
    IntervalNode *y = x->right;
    IntervalNode *T2 = y->left;

    // Perform rotation
    y->left = x;
    x->right = T2;

    // Update heights and max_end
    updateHeight(x);
    updateHeight(y);
    updateMaxEnd(x);
    updateMaxEnd(y);

    return y;
  }

  IntervalNode *getParent(IntervalNode *node) {
    if (!node || node == root)
      return nullptr;
    return findParent(root, node);
  }

  IntervalNode *findParent(IntervalNode *current, IntervalNode *target) {
    if (!current || current == target)
      return nullptr;
    if (current->left == target || current->right == target)
      return current;

    IntervalNode *leftResult = findParent(current->left, target);
    if (leftResult)
      return leftResult;
    return findParent(current->right, target);
  }

  void destroyTree(IntervalNode *node) {
    if (node) {
      destroyTree(node->left);
      destroyTree(node->right);
      delete node;
    }
  }

  void inorderTraversal(IntervalNode *node,
                        std::vector<IntervalNode *> &result) {
    if (node) {
      inorderTraversal(node->left, result);
      result.push_back(node);
      inorderTraversal(node->right, result);
    }
  }

  IntervalNode *insertNodeBalanced(IntervalNode *node, IntervalNode *newNode) {
    if (!node)
      return newNode;

    // Insert based on start position (and end as tiebreaker)
    if (newNode->start < node->start ||
        (newNode->start == node->start && newNode->end < node->end)) {
      node->left = insertNodeBalanced(node->left, newNode);
    } else {
      node->right = insertNodeBalanced(node->right, newNode);
    }

    // Update height and max_end
    updateHeight(node);
    updateMaxEnd(node);

    // Get balance factor
    int balance = getBalance(node);

    // Case 1: Left Left
    if (balance > 1 && (newNode->start < node->left->start ||
                        (newNode->start == node->left->start &&
                         newNode->end < node->left->end))) {
      return rotateRight(node);
    }

    // Case 2: Right Right
    if (balance < -1 && (newNode->start > node->right->start ||
                         (newNode->start == node->right->start &&
                          newNode->end > node->right->end))) {
      return rotateLeft(node);
    }

    // Case 3: Left Right
    if (balance > 1 && (newNode->start > node->left->start ||
                        (newNode->start == node->left->start &&
                         newNode->end > node->left->end))) {
      node->left = rotateLeft(node->left);
      return rotateRight(node);
    }

    // Case 4: Right Left
    if (balance < -1 && (newNode->start < node->right->start ||
                         (newNode->start == node->right->start &&
                          newNode->end < node->right->end))) {
      node->right = rotateRight(node->right);
      return rotateLeft(node);
    }

    return node;
  }

public:
  IntervalTree() : root(nullptr) {}

  ~IntervalTree() { destroyTree(root); }

  void insert(uint32_t start, uint32_t end,
              const tid_t &tid) {
    
    //GMessage("start = %d\t end = %d\t tid = %d\n", start, end, tid);
    // Insert tid to set of all TIDs
    all_tids.insert(tid);

    // Find tree nodes that overlap
    auto overlapping = findAllOverlapping(start, end);

    // No overlaps: create new node for entire range
    if (overlapping.empty()) {
      auto newNode = new IntervalNode(start, end);
      newNode->tids.insert(tid);
      insertNode(newNode);
      return;
    }

    // Check for exact match
    for (auto &node : overlapping) {
      if (node->start == start && node->end == end) {
        node->tids.insert(tid);
        return;
      }
    }

    // *^*^*^ *^*^*^ *^*^*^ *^*^*^
    // Split nodes
    // *^*^*^ *^*^*^ *^*^*^ *^*^*^

    std::set<uint32_t> split_points;
    for (auto &node : overlapping) {
      if (node->start < start && start < node->end) {
        split_points.insert(start);
      }
      if (node->start < end && end < node->end) {
        split_points.insert(end);
      }
    }

    // Perform all splits in one pass
    for (uint32_t split_point : split_points) {
      auto nodes_at_point = findAllOverlapping(split_point, split_point + 1);
      for (auto &node : nodes_at_point) {
        if (node->start < split_point && node->end > split_point) {
          splitInterval(node, split_point);
        }
      }
    }

    // Now find all intervals that should contain this TID
    auto final_overlapping = findAllOverlapping(start, end);

    // Track which parts of [start, end) are covered by existing nodes
    std::vector<std::pair<uint32_t, uint32_t>> covered_ranges;
    covered_ranges.reserve(final_overlapping.size());

    for (auto &node : final_overlapping) {
      uint32_t overlap_start = std::max(node->start, start);
      uint32_t overlap_end = std::min(node->end, end);

      if (overlap_start < overlap_end) {
        node->tids.insert(tid);
        covered_ranges.emplace_back(overlap_start, overlap_end);
      }
    }

    // *^*^*^ *^*^*^ *^*^*^ *^*^*^
    // Create nodes for gaps
    // *^*^*^ *^*^*^ *^*^*^ *^*^*^

    std::sort(covered_ranges.begin(), covered_ranges.end());
    auto merged_ranges = mergeRanges(covered_ranges);
    std::vector<IntervalNode *> new_nodes;

    // There's a gap at beginning
    if (!merged_ranges.empty() && merged_ranges[0].first > start) {
      new_nodes.push_back(new IntervalNode(start, merged_ranges[0].first));
    }

    // There's a gap between ranges
    for (size_t i = 0; i + 1 < merged_ranges.size(); i++) {
      if (merged_ranges[i].second < merged_ranges[i + 1].first) {
        new_nodes.push_back(new IntervalNode(merged_ranges[i].second,
                                            merged_ranges[i + 1].first));
      }
    }

    // There's a gap at end
    if (!merged_ranges.empty() && merged_ranges.back().second < end) {
      new_nodes.push_back(new IntervalNode(merged_ranges.back().second, end));
    }

    // Handle case where no ranges were covered
    if (merged_ranges.empty()) {
      new_nodes.push_back(new IntervalNode(start, end));
    }

    // Batch insert new nodes
    for (auto *node : new_nodes) {
      node->tids.insert(tid);
      insertNode(node);
    }
  }

private:
  void insertNode(IntervalNode *newNode) {
    root = insertNodeBalanced(root, newNode);
  }

  std::vector<IntervalNode *> findAllOverlapping(uint32_t start, uint32_t end) {
    std::vector<IntervalNode *> result;
    findOverlappingHelper(root, start, end, result);
    return result;
  }

  void findOverlappingHelper(IntervalNode *node, uint32_t start, uint32_t end,
                             std::vector<IntervalNode *> &result) {
    if (!node)
      return;

    // Check if current interval overlaps
    if (node->start < end && node->end > start) {
      result.push_back(node);
    }

    // Search left subtree
    if (node->left && node->left->max_end > start) {
      findOverlappingHelper(node->left, start, end, result);
    }

    // Search right subtree
    if (node->right && node->start < end) {
      findOverlappingHelper(node->right, start, end, result);
    }
  }

  // Merge several ranges together
  std::vector<std::pair<uint32_t, uint32_t>>
  mergeRanges(std::vector<std::pair<uint32_t, uint32_t>> &ranges) {
    if (ranges.empty())
      return ranges;

    std::vector<std::pair<uint32_t, uint32_t>> merged;
    merged.reserve(ranges.size());
    merged.push_back(ranges[0]);

    for (size_t i = 1; i < ranges.size(); i++) {
      // Overlapping or adjacent
      if (ranges[i].first <= merged.back().second) {
        merged.back().second = std::max(merged.back().second, ranges[i].second);
      } else {
        merged.push_back(ranges[i]);
      }
    }

    return merged;
  }

  // Split interval in two
  void splitInterval(IntervalNode *node, uint32_t split_point) {
    std::set<tid_t> tids = node->tids; // std::string> tids = node->tids;

    // Create new node for second half
    auto newNode = new IntervalNode(split_point, node->end);
    newNode->tids = tids;

    node->end = split_point; // update original node to represent first half
    insertNode(newNode);     // insert the new node into tree
  }

  void getTidNodesHelper(IntervalNode *node,
                         const tid_t &tid, // const std::string &tid,
                         std::vector<IntervalNode *> &result) {
    if (!node)
      return;

    if (node->tids.count(tid)) {
      result.push_back(node);
    }

    getTidNodesHelper(node->left, tid, result);
    getTidNodesHelper(node->right, tid, result);
  }

  // Get all nodes containing a specific TID
  std::vector<IntervalNode *>
  getTidNodes(const tid_t &tid) { // std::string &tid) {
    std::vector<IntervalNode *> result;
    getTidNodesHelper(root, tid, result);
    return result;
  }

  // Build chain for a specific TID
  void buildTidChain(const tid_t &tid) { // const std::string &tid) {
    // Get all nodes with this TID in genomic order
    auto tidNodes = getTidNodes(tid);

    // Note: around 15s faster to sort here than use tid_chains

    // Sort by genomic position
    std::sort(tidNodes.begin(), tidNodes.end(),
              [](IntervalNode *a, IntervalNode *b) {
                if (a->start != b->start)
                  return a->start < b->start;
                return a->end < b->end;
              });

    // Add to tid_chains
    //tid_chains[tid] = tidNodes;

    // Build bidirectional chain
    for (size_t i = 0; i < tidNodes.size(); ++i) {
      if (i > 0) {
        tidNodes[i]->tid_prev_nodes[tid] = tidNodes[i - 1];
        tidNodes[i - 1]->tid_next_nodes[tid] = tidNodes[i];
      }
    }
  }

  // Precompute cumulative lengths so that match positions can be calculated
  void precomputeCumulativeLengths(const tid_t &tid) { // std::string &tid) {
    //auto tidChain = tid_chains[tid];
    auto tidChain = getTidNodes(tid);

    // Need to sort by genomic position
    // TODO: store these chains so they don't need to be resorted
    std::sort(tidChain.begin(), tidChain.end(),
              [](IntervalNode *a, IntervalNode *b) {
                if (a->start != b->start)
                  return a->start < b->start;
                return a->end < b->end;
              });

    uint32_t cumulative = 0;
    uint32_t prev_end = -1;
    for (size_t i = 0; i < tidChain.size(); ++i) {
      auto *node = tidChain[i];

      // Subtract 1 if this and previous node are directly adjacent
      if (node->start == prev_end) {
        node->tid_cum_len[tid] = cumulative - 1;
        cumulative += (node->end - node->start);
      } else {
        node->tid_cum_len[tid] = cumulative;
        cumulative += (node->end - node->start) + 1;
      }

      prev_end = node->end;
    }
  }

public:
  // For debugging (print_tree)
  std::vector<IntervalNode *> getOrderedIntervals() {
    std::vector<IntervalNode *> result;
    inorderTraversal(root, result);
    return result;
  }

  // Find all intervals that overlap with the given range
  std::vector<IntervalNode *> findOverlapping(uint32_t start, uint32_t end) {
    std::vector<IntervalNode *> result;
    findOverlappingHelper(root, start, end, result);
    std::sort(
        result.begin(), result.end(),
        [](IntervalNode *a, IntervalNode *b) { return a->start < b->start; });
    return result;
  }

  // Get next node in chain for a specific TID
  IntervalNode *
  getNextNodeForTid(IntervalNode *node,
                    const tid_t &tid) { // const std::string &tid) {
    // if (tid_chains.find(tid) == tid_chains.end()) {
    //   buildTidChain(tid);
    // }
    auto it = node->tid_next_nodes.find(tid);
    return (it != node->tid_next_nodes.end()) ? it->second : nullptr;
  }

  // Get previous node in chain for a specific TID
  IntervalNode *
  getPrevNodeForTid(IntervalNode *node,
                    const tid_t &tid) { // const std::string &tid) {
    // if (tid_chains.find(tid) == tid_chains.end()) {
    //   buildTidChain(tid);
    // }
    auto it = node->tid_prev_nodes.find(tid);
    return (it != node->tid_prev_nodes.end()) ? it->second : nullptr;
  }

  // Build all TID chains after tree construction is complete
  void buildAllTidChains() {
    for (const auto &tid : all_tids) {
      buildTidChain(tid);
    }
  }

  // Precompute for all TIDs in the tree
  void precomputeAllCumulativeLengths() {
    for (const auto &tid : all_tids) {
      precomputeCumulativeLengths(tid);
    }
  }

  // Check for cumulative length
  uint32_t findCumulativeLength(IntervalNode *node, const tid_t &tid) {
    auto it = node->tid_cum_len.find(tid);
    return (it != node->tid_cum_len.end()) ? it->second : 0;
  }
};

// g2t tree using interval tree
struct g2tTree {
  // map from transcript names to transcript IDs
  std::unordered_map<tid_t, std::string> name_id_map;
  std::vector<std::string> tid_names;
  const std::string invalid_name{"INVALID TID"};

  IntervalTree *fw_tree; // tree for forward strand guides
  IntervalTree *rc_tree; // tree for reverse strand guides

  g2tTree() {
    fw_tree = new IntervalTree();
    rc_tree = new IntervalTree();
  }

  ~g2tTree() {
    delete fw_tree;
    delete rc_tree;
  }

private:
  IntervalTree *getTreeForStrand(char strand) {
    if (strand == '+' || strand == 1)
      return fw_tree;
    else if (strand == '-' || strand == -1)
      return rc_tree;
    else
      return nullptr;
  }

public:
  // If tid_name is already part of the tree, then return the existing
  // tid for it, otherwise give it the next available tid and return that.
  tid_t insertTidString(const char *&tid_name, BamIO* io) {
    //tid_t next_tid = static_cast<tid_t>(name_id_map.size());
    tid_t tid = io->get_tid(tid_name);

    // auto lookup = name_id_map.find(tid_name);
    // if (lookup == name_id_map.end()) {
    //   name_id_map.insert({tid_name, next_tid});
    //   tid_names.push_back(tid_name);
    //   return next_tid;

    // } else {
    //   return lookup->second;
    // }

    auto lookup = name_id_map.find(tid);
    if (lookup == name_id_map.end()) {
      name_id_map.insert({tid, tid_name});
      tid_names.push_back(tid_name);
    }

    return tid;
  }

  // get the string name given an id
  const std::string &getTidName(tid_t id) {
    return name_id_map.find(id) != name_id_map.end() ? name_id_map[id] : invalid_name;
  }

  // Add guide exon with TID and transcript start
  void addGuideExon(uint32_t start, uint32_t end,
                    const tid_t &tid, // const std::string &tid,
                    char strand) {
    IntervalTree *tree = getTreeForStrand(strand);
    if (tree)
      tree->insert(start, end, tid);
  }

  // Call this after loading all transcript data
  void buildAllTidChains() {
    fw_tree->buildAllTidChains();
    rc_tree->buildAllTidChains();
  }

  // Call this after loading all transcript data
  // TODO: only build when first see TID in read
  void precomputeAllCumulativeLengths() {
    fw_tree->precomputeAllCumulativeLengths();
    rc_tree->precomputeAllCumulativeLengths();
  }

  // Find all guide TIDs that overlap with a read exon
  std::vector<IntervalNode *> getIntervals(uint32_t readStart, uint32_t readEnd,
                                           char strand) {
    IntervalTree *tree = getTreeForStrand(strand);
    if (!tree)
      return std::vector<IntervalNode *>();

    return tree->findOverlapping(readStart, readEnd);
  }

  // Get cumulative previous size of exons from transcript
  uint32_t getCumulativeLength(IntervalNode *node, const tid_t &tid,
                               char strand) {
    IntervalTree *tree = getTreeForStrand(strand);
    return tree->findCumulativeLength(node, tid);
  }

  // Get next node in chain for a specific TID
  IntervalNode *getNextNode(IntervalNode *node,
                            const tid_t &tid, //) {const std::string &tid,
                            char strand) {
    IntervalTree *tree = getTreeForStrand(strand);
    return tree->getNextNodeForTid(node, tid);
  }

  // Get previous node in chain for a specific TID
  IntervalNode *getPrevNode(IntervalNode *node,
                            const tid_t &tid, // const std::string &tid,
                            char strand) {
    IntervalTree *tree = getTreeForStrand(strand);
    return tree->getPrevNodeForTid(node, tid);
  }
};

struct ReadInfo {
  std::set<std::tuple<tid_t, uint32_t>> matches;
  GSamRecord *brec;
  bool valid_read;
  uint32_t read_index;
  uint32_t read_size;
  uint32_t nh_i;
   bool is_paired;
  bool has_introns;
  bool is_reverse;
  uint32_t soft_clip_front;
  uint32_t soft_clip_back;

  ReadInfo() : matches(), brec(nullptr), valid_read(false), read_index(0), 
               read_size(0), nh_i(0), is_paired(false), has_introns(false), 
               is_reverse(false), soft_clip_front(0), soft_clip_back(0) {}

  ~ReadInfo() {
    // brec is deleted by CReadAln
  }
};

struct BamInfo {
  bool same_transcript; // true if both mates map to same transcript
  bool valid_pair;      // true if mate information should be output
  bool is_paired;       // true if from paired reads
  
  // Read 1 information
  read_id_t read_index;
  tid_t tid;
  uint32_t pos;
  uint32_t nh;
  uint32_t read_size;
  GSamRecord *brec;
  bool has_introns;
  bool is_reverse;
  uint32_t soft_clip_front;
  uint32_t soft_clip_back;
  bool discard_read;

  // Read 2 information
  read_id_t mate_index;
  tid_t mate_tid;
  uint32_t mate_pos;
  uint32_t mate_nh;
  uint32_t mate_size;
  GSamRecord *mate_brec;
  bool mate_has_introns;
  bool mate_is_reverse;
  uint32_t mate_soft_clip_front;
  uint32_t mate_soft_clip_back;
  bool mate_discard_read;

  BamInfo() : discard_read(false), mate_discard_read(false) {}

  ~BamInfo() {
    // idk who deletes this but its somebody
  }
};

// Collect all refguide transcripts for a single genomic sequence
struct GRefData {
  GList<GffObj> rnas; // all transcripts on this genomic seq
  int gseq_id;
  const char *gseq_name;
  GRefData(int gid = -1)
      : rnas(false, false, false), gseq_id(gid), gseq_name(NULL) {
    gseq_id = gid;
    if (gseq_id >= 0)
      gseq_name = GffObj::names->gseqs.getName(gseq_id);
  }

  void add(GffReader *gffr, GffObj *t) {
    if (gseq_id >= 0) {
      if (gseq_id != t->gseq_id)
        GError("Error: invalid call to GRefData::add() - different genomic "
               "sequence!\n");
    } else { // adding first transcript, initialize storage
      gseq_id = t->gseq_id;
      gseq_name = t->getGSeqName();
      if (gffr->gseqtable[gseq_id] == NULL)
        GError("Error: invalid genomic sequence data (%s)!\n", gseq_name);
      rnas.setCapacity(gffr->gseqtable[gseq_id]->fcount);
    }
    rnas.Add(t);
    t->isUsed(true);
    // setLocus(t); //use the GRefLocus::mexons to quickly find an overlap with
    // existing loci, or create a new one
  }

  bool operator==(GRefData &d) { return gseq_id == d.gseq_id; }
  bool operator<(GRefData &d) { return (gseq_id < d.gseq_id); }
};

// Bundle status
enum BundleStatus {
  BUNDLE_STATUS_CLEAR = 0, // Available for loading/prepping
  BUNDLE_STATUS_LOADING = 1, // Being prepared by the main thread (there can be only one)
  BUNDLE_STATUS_READY = 2 // Ready to be processed, or being processed
};

// Read Alignment
struct CReadAln : public GSeg {
  char strand; // 1, 0 (unknown), -1 (reverse)
  short int nh;
  uint32_t len;
  float read_count; // keeps count for all reads (including paired and unpaired)
  bool unitig : 1;  // set if read come from an unitig
  bool longread : 1;      // set if read comes from long read data
  GVec<float> pair_count; // keeps count for all paired reads
  GVec<int> pair_idx; // keeps indeces for all pairs in assembly mode, or all
                      // reads that were collapsed in merge mode

  // For mate position calculations
  uint32_t mate_genomic_pos;
  tid_t mate_transcript_id;
  bool mate_found;

  GVec<GSeg> segs; //"exons"
  GSamRecord *brec; // store BAM record

  CReadAln(char _strand = 0, short int _nh = 0, int rstart = 0, int rend = 0)
      : GSeg(rstart, rend), // name(rname),
        strand(_strand), nh(_nh), len(0), read_count(0), unitig(false),
        longread(false), pair_count(), pair_idx(), mate_genomic_pos(0),
        mate_transcript_id(), mate_found(false), segs(), brec() {}

  CReadAln(CReadAln &rd) : GSeg(rd.start, rd.end) { // copy contructor
    strand = rd.strand;
    nh = rd.nh;
    len = rd.len;
    read_count = rd.read_count;
    unitig = rd.unitig;
    longread = rd.longread;
    pair_count = rd.pair_count;
    pair_idx = rd.pair_idx;
  }

  ~CReadAln() { delete brec; }

  friend std::ostream &operator<<(std::ostream &, const CReadAln &);
};

inline std::ostream &operator<<(std::ostream &os, const CReadAln &aln) {
  os << "[CReadAln](strand: " << aln.strand << ", nh: " << aln.nh
     << ", len: " << aln.len << ", read_count: " << aln.read_count
     << ", unitig: " << aln.unitig << ", long_read: " << aln.longread
     << ", mate_genomic_pos: " << aln.mate_genomic_pos
     << ", mate_transcript_id: " << aln.mate_transcript_id
     << ", mate_found: " << aln.mate_found << ")";
  return os;
}

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

  AlnGroups() : key_to_group(), groups() {}

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

// Bundle data structure, holds input data parsed from BAM file
struct BundleData {
  BundleStatus status; // bundle status

  int idx; // index in the main bundles array
  int start; // bundle start
  int end; // bundle end
  unsigned long num_reads; // number of reads in this bundle
  double num_fragments;    // aligned read/pairs

  GStr refseq;           // reference sequence name
  char *gseq;            // genomic sequence for the bundle
  GList<CReadAln> reads; // all reads in bundle
  GPVec<GffObj> guides; // all guides in bundle

BundleData():status(BundleStatus::BUNDLE_STATUS_CLEAR), idx(0), start(0), end(0), 
			refseq(), gseq(NULL), reads(false, true), guides(false) //rc_data(NULL) { } // sorted, free elements
			{}
 
  // Called when the bundle is valid and ready to be processed
  void getReady(int curr_start, int curr_end) {
    start = curr_start;
    end = curr_end;
    status = BundleStatus::BUNDLE_STATUS_READY;
  }

  void keepGuide(GffObj *t) {
    guides.Add(t);
  }

  void Clear() {

    guides.Clear();
    reads.Clear();
    // TODO: Think about managing this with a unique_ptr
    GFREE(gseq);

    start = 0;
    end = 0;
    status = BundleStatus::BUNDLE_STATUS_CLEAR;
  }

  ~BundleData() { Clear(); }
};

struct WorkerArgs {
  GPVec<BundleData> *bundle_queue;
  BamIO *io;

  WorkerArgs(GPVec<BundleData> *q, BamIO *d) : bundle_queue(q), io(d) {}
};

#endif
