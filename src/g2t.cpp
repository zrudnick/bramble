
#include <numeric>
#include <algorithm>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "types.h"
#include "bundles.h"
#include "bramble.h"
#include "g2t.h"


namespace bramble {

  IntervalNode::IntervalNode(uint32_t s, uint32_t e)
    : start(s), 
      end(e), 
      max_end(e), 
      height(1), 
      left(nullptr), 
      right(nullptr) {}

  IntervalTree::IntervalTree() : root(nullptr) {}

  IntervalTree::~IntervalTree() { destroyTree(root); }
  
  int IntervalTree::getHeight(IntervalNode *node) { 
    return node ? node->height : 0; 
  }

  int IntervalTree::getBalance(IntervalNode *node) {
    return node ? getHeight(node->left) - getHeight(node->right) : 0;
  }

  void IntervalTree::updateHeight(IntervalNode *node) {
    if (node) {
      uint32_t left_height = getHeight(node->left);
      uint32_t right_height = getHeight(node->right);
      node->height = 1 + std::max(left_height, right_height);
    }
  }

  void IntervalTree::updateMaxEnd(IntervalNode *node) {
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

  IntervalNode *IntervalTree::rotateRight(IntervalNode *y) {

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

  IntervalNode *IntervalTree::rotateLeft(IntervalNode *x) {
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

  IntervalNode *IntervalTree::getParent(IntervalNode *node) {
    if (!node || node == root)
      return nullptr;
    return findParent(root, node);
  }

  IntervalNode 
  *IntervalTree::findParent(IntervalNode *current, IntervalNode *target) {
    if (!current || current == target)
      return nullptr;
    if (current->left == target || current->right == target)
      return current;

    IntervalNode *leftResult = findParent(current->left, target);
    if (leftResult)
      return leftResult;
    return findParent(current->right, target);
  }

  void IntervalTree::destroyTree(IntervalNode *node) {
    if (node) {
      destroyTree(node->left);
      destroyTree(node->right);
      delete node;
    }
  }

  void IntervalTree::inorderTraversal(IntervalNode *node,
                                      std::vector<IntervalNode *> &result) {
    if (node) {
      inorderTraversal(node->left, result);
      result.push_back(node);
      inorderTraversal(node->right, result);
    }
  }

  IntervalNode 
  *IntervalTree::insertNodeBalanced(IntervalNode *node, 
                                    IntervalNode *newNode) {
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

  void IntervalTree::insert(uint32_t start, uint32_t end, 
                            const tid_t &tid) {

    if (start == end) return;
    
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

      if (node->start >= start && node->end <= end) {
        node->tids.insert(tid);
        covered_ranges.emplace_back(node->start, node->end);
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

  void IntervalTree::insertNode(IntervalNode *newNode) {
    root = insertNodeBalanced(root, newNode);
  }

  std::vector<IntervalNode *> 
  IntervalTree::findAllOverlapping(uint32_t start, uint32_t end) {
    std::vector<IntervalNode *> result;
    findOverlappingHelper(root, start, end, result);
    return result;
  }

  void IntervalTree::findOverlappingHelper(IntervalNode *node, 
                                          uint32_t start, uint32_t end,
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
  IntervalTree::mergeRanges(std::vector<std::pair<uint32_t, uint32_t>> &ranges) {
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
  void IntervalTree::splitInterval(IntervalNode *node, 
                                  uint32_t split_point) {
    std::set<tid_t> tids = node->tids;

    // Create new node for second half
    auto newNode = new IntervalNode(split_point, node->end);
    newNode->tids = tids;

    node->end = split_point; // update original node to represent first half
    insertNode(newNode);     // insert the new node into tree
  }

  void IntervalTree::getTidNodesHelper(IntervalNode *node,
                                      const tid_t &tid,
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
  IntervalTree::getTidNodes(const tid_t &tid) {
    std::vector<IntervalNode *> result;
    getTidNodesHelper(root, tid, result);
    return result;
  }

  // Build chain for a specific TID
  void IntervalTree::buildTidChain(const tid_t &tid) {
    // Get all nodes with this TID in genomic order
    auto tidNodes = getTidNodes(tid);

    // Sort by genomic position
    std::sort(tidNodes.begin(), tidNodes.end(),
              [](IntervalNode *a, IntervalNode *b) {
                if (a->start != b->start)
                  return a->start < b->start;
                return a->end < b->end;
              });

    // Build bidirectional chain
    for (size_t i = 0; i < tidNodes.size(); ++i) {
      if (i > 0) {
        tidNodes[i]->tid_prev_nodes[tid] = tidNodes[i - 1];
        tidNodes[i - 1]->tid_next_nodes[tid] = tidNodes[i];
      }
    }
  }

  void IntervalTree::precomputeCumulativeLengths(const tid_t &tid, char strand) {
    auto tidChain = getTidNodes(tid);

    // Need to sort by genomic position
    std::sort(tidChain.begin(), tidChain.end(),
              [](IntervalNode *a, IntervalNode *b) {
                if (a->start != b->start)
                  return a->start < b->start;
                return a->end < b->end;
              });

    uint32_t total_length = 0;
    for (size_t i = 0; i < tidChain.size(); ++i) {
      auto *node = tidChain[i];
      total_length += (node->end - node->start);
    }

    uint32_t cumulative = 0;
    for (size_t i = 0; i < tidChain.size(); ++i) {
      auto *node = tidChain[i];

      uint32_t node_length = node->end - node->start;
      uint32_t left = cumulative;
      uint32_t right =  total_length - cumulative - node_length;
      node->tid_cum_len[tid] = left;
      cumulative += node_length;
    }
    transcript_lengths[tid] = cumulative;
  }

  // For debugging (print_tree)
  std::vector<IntervalNode *> IntervalTree::getOrderedIntervals() {
    std::vector<IntervalNode *> result;
    inorderTraversal(root, result);
    return result;
  }

  // Find all intervals that overlap with the given range
  std::vector<IntervalNode *> 
  IntervalTree::findOverlapping(uint32_t start, uint32_t end) {
    std::vector<IntervalNode *> result;
    findOverlappingHelper(root, start, end, result);
    std::sort(
        result.begin(), result.end(),
        [](IntervalNode *a, IntervalNode *b) { return a->start < b->start; });

    return result;
  }

  // Get next node in chain for a specific TID
  IntervalNode *
  IntervalTree::getNextNodeForTid(IntervalNode *node, const tid_t &tid) {
    auto it = node->tid_next_nodes.find(tid);
    return (it != node->tid_next_nodes.end()) ? it->second : nullptr;
  }

  // Get previous node in chain for a specific TID
  IntervalNode *
  IntervalTree::getPrevNodeForTid(IntervalNode *node, const tid_t &tid) {
    auto it = node->tid_prev_nodes.find(tid);
    return (it != node->tid_prev_nodes.end()) ? it->second : nullptr;
  }

  // Build all TID chains after tree construction is complete
  void IntervalTree::buildAllTidChains() {
    for (const auto &tid : all_tids) {
      buildTidChain(tid);
    }
  }

  // Precompute for all TIDs in the tree
  void IntervalTree::precomputeAllCumulativeLengths(char strand) {
    for (const auto &tid : all_tids) {
      precomputeCumulativeLengths(tid, strand);
    }
  }

  // Check for cumulative length
  uint32_t 
  IntervalTree::findCumulativeLength(IntervalNode *node, 
                                    const tid_t &tid) {
    auto it = node->tid_cum_len.find(tid);
    return (it != node->tid_cum_len.end()) ? it->second : 0;
  }

  uint32_t IntervalTree::getTranscriptLength(tid_t tid) {
    return transcript_lengths[tid];
  }

  g2tTree::g2tTree() {
    fw_tree = new IntervalTree();
    rc_tree = new IntervalTree();
  }

  g2tTree::~g2tTree() {
    delete fw_tree;
    delete rc_tree;
  }

  IntervalTree *g2tTree::getTreeForStrand(char strand) {
    if (strand == '+' || strand == 1)
      return fw_tree;
    else if (strand == '-' || strand == -1)
      return rc_tree;
    else
      return nullptr;
  }

  // If tid_name is already part of the tree, then return the existing
  // tid for it, otherwise give it the next available tid and return that.
  tid_t g2tTree::insertTidString(const char *&tid_name, BamIO* io) {
    tid_t tid = io->get_tid(tid_name);

    auto lookup = name_id_map.find(tid);
    if (lookup == name_id_map.end()) {
      name_id_map.insert({tid, tid_name});
      tid_names.push_back(tid_name);
    }

    return tid;
  }

  // get the string name given an id
  const std::string &g2tTree::getTidName(tid_t id) {
    return (name_id_map.find(id) != name_id_map.end()) ? 
      name_id_map[id] : invalid_name;
  }

  // Add guide exon with TID and transcript start
  void g2tTree::addGuideExon(uint32_t start, uint32_t end,
                            const tid_t &tid, char strand) {
    IntervalTree *tree = getTreeForStrand(strand);
    if (tree)
      tree->insert(start, end, tid);
  }

  // Call this after loading all transcript data
  void g2tTree::buildAllTidChains() {
    fw_tree->buildAllTidChains();
    rc_tree->buildAllTidChains();
  }

  // Call this after loading all transcript data
  void g2tTree::precomputeAllCumulativeLengths() {
    fw_tree->precomputeAllCumulativeLengths('+');
    rc_tree->precomputeAllCumulativeLengths('-');
  }

  // Find all guide TIDs that overlap with a read exon
  std::vector<IntervalNode *> 
  g2tTree::getIntervals(uint32_t readStart, uint32_t readEnd,
                        char strand) {
    IntervalTree *tree = getTreeForStrand(strand);
    if (!tree)
      return std::vector<IntervalNode *>();

    return tree->findOverlapping(readStart, readEnd);
  }

  // Get cumulative previous size of exons from transcript
  uint32_t 
  g2tTree::getCumulativeLength(IntervalNode *node, 
                              const tid_t &tid, char strand) {
    IntervalTree *tree = getTreeForStrand(strand);
    return tree->findCumulativeLength(node, tid);
  }

  // Get next node in chain for a specific TID
  IntervalNode *g2tTree::getNextNode(IntervalNode *node, const tid_t &tid, 
                                    char strand) {
    IntervalTree *tree = getTreeForStrand(strand);
    return tree->getNextNodeForTid(node, tid);
  }

  // Get previous node in chain for a specific TID
  IntervalNode *g2tTree::getPrevNode(IntervalNode *node, const tid_t &tid,
                                    char strand) {
    IntervalTree *tree = getTreeForStrand(strand);
    return tree->getPrevNodeForTid(node, tid);
  }

  // Get previous node in chain for a specific TID
  uint32_t g2tTree::getTranscriptLength(const tid_t &tid,
                                            char strand) {
    IntervalTree *tree = getTreeForStrand(strand);
    return tree->getTranscriptLength(tid);
  }
  
};