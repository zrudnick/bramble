
#pragma once
#include <algorithm>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "GArgs.h"
#include "GBitVec.h"
#include "GHashMap.hh"
#include "GSam.h"
#include "GStr.h"
#include "gff.h"
#include "time.h"
#include "types.h"

using namespace bramble;

namespace bramble {

struct BamIO;

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

  void printTreeStructure() {
    std::cout << "\n=== TREE STRUCTURE ===" << std::endl;
    if (!root) {
        std::cout << "Tree is empty" << std::endl;
        return;
    }
    printNodeRecursive(root, "", true, 0);
    std::cout << "=====================\n" << std::endl;
}

void printNodeRecursive(IntervalNode *node, const std::string &prefix, bool isLast, int depth) {
    if (!node) {
        std::cout << prefix << (isLast ? "└── " : "├── ") << "NULL" << std::endl;
        return;
    }
    
    // Print current node
    std::cout << prefix << (isLast ? "└── " : "├── ") 
              << "[" << node->start << "," << node->end << ") "
              << "h=" << node->height 
              << " max=" << node->max_end
              << " tids={";
    
    // Print first few TIDs
    int count = 0;
    for (auto tid : node->tids) {
        if (count++ > 0) std::cout << ",";
        std::cout << tid;
        if (count >= 3) { std::cout << "..."; break; }
    }
    std::cout << "}" << std::endl;
    
    // Print children if they exist OR if we want to show NULL children for debugging
    bool hasChildren = (node->left != nullptr || node->right != nullptr);
    
    if (hasChildren || depth < 3) {  // Show NULL children for shallow depths
        std::string newPrefix = prefix + (isLast ? "    " : "│   ");
        
        if (node->left || node->right) {
            printNodeRecursive(node->left, newPrefix, !node->right, depth + 1);
            if (node->right) {
                printNodeRecursive(node->right, newPrefix, true, depth + 1);
            }
        }
    }
  }

  IntervalNode *rotateRight(IntervalNode *y) {
    if (!y->left) {
        std::cout << "ERROR: rotateRight called but y->left is NULL!" << std::endl;
        std::cout << "Node y: [" << y->start << "," << y->end << "] has no left child" << std::endl;
        std::cout << "Tree structure before failed rotation:" << std::endl;
        printTreeStructure();
        return y;  // Return unchanged rather than segfault
    }

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

    if (start == end) return;
    
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

}