#ifndef MAIN_H
#define MAIN_H

#include "GArgs.h"
#include "GStr.h"
#include "gff.h"
#include "GSam.h"
#include "GBitVec.h"
#include "time.h"
#include "GHashMap.hh"
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <memory>
#include <iostream>
using namespace std;

#define MAX_NODE 1000000
#define SMALL_EXON 35 			// exons smaller than this have a tendency to be missed by long read data

#define GMEMTRACE 0				// Memory tracing
#ifndef NOTHREADS				// Use threads
#define THREADING_ENABLED 1
#else
#define THREADING_ENABLED 0
#endif

#define RC_MIN_EOVL 5

constexpr uint32_t longintron =
     100000; // don't trust introns longer than this unless there is higher
             // evidence; about 99% of all human annotated introns are shorter
             // than this
constexpr uint32_t longintronanchor = 25;
constexpr float mismatchfrac = 0.02;
using tid_t = uint32_t;

extern bool verbose;
extern bool debugMode;

struct IntervalNode {
   uint32_t start, end;
   uint32_t max_end; // Maximum end value in subtree (interval tree augmentation)
   uint32_t height;
   std::set<tid_t> tids;
   // does this need to be an ordered map?
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
 
  void insert(uint32_t start, uint32_t end, const tid_t& tid) { //const std::string &tid) {
     // Insert tid to set of all TIDs
     all_tids.insert(tid);
 
     // Find tree nodes that overlap
     auto overlapping = findAllOverlapping(start, end);
 
     // No overlaps: create new node for entire range
     if (overlapping.empty()) {
       auto newNode = new IntervalNode(start, end);
       newNode->tids.insert(tid);
       // newNode->height = 1;
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
     for (size_t i = 0; i < merged_ranges.size() - 1; i++) {
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
       // node->height = 1;
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
    std::set<tid_t> tids = node->tids;//std::string> tids = node->tids;
 
     // Create new node for second half
     auto newNode = new IntervalNode(split_point, node->end);
     newNode->tids = tids;
     // newNode->height = 1;
 
     node->end = split_point; // update original node to represent first half
     insertNode(newNode);     // insert the new node into tree
   }
 
  void getTidNodesHelper(IntervalNode *node, const tid_t &tid, //const std::string &tid,
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
  std::vector<IntervalNode *> getTidNodes(const tid_t &tid) { //std::string &tid) {
     std::vector<IntervalNode *> result;
     getTidNodesHelper(root, tid, result);
     return result;
   }
 
   // Build chain for a specific TID
  void buildTidChain(const tid_t &tid) { //const std::string &tid) {
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
  IntervalNode *getNextNodeForTid(IntervalNode *node, const tid_t& tid) { //const std::string &tid) {
     auto it = node->tid_next_nodes.find(tid);
     return (it != node->tid_next_nodes.end()) ? it->second : nullptr;
   }
 
   // Get previous node in chain for a specific TID
  IntervalNode *getPrevNodeForTid(IntervalNode *node, const tid_t& tid) { //const std::string &tid) {
     auto it = node->tid_prev_nodes.find(tid);
     return (it != node->tid_prev_nodes.end()) ? it->second : nullptr;
   }
 
   // Build all TID chains after tree construction is complete
   void buildAllTidChains() {
     for (const auto &tid : all_tids) {
       buildTidChain(tid);
     }
   }
 
   // Precompute cumulative lengths so that match positions can be calculated
  void precomputeCumulativeLengths(const tid_t &tid) { //std::string &tid) {
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
  std::unordered_map<std::string, tid_t> name_id_map;
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

    IntervalTree* getTreeForStrand(char strand) {
        if (strand == '+' || strand == 1) return fw_tree;
        else if (strand == '-' || strand == -1) return rc_tree;
        else return nullptr;
    }
    
public:

  // if tid_name is already part of the tree, then return the existing 
  // tid for it, otherwise give it the next available tid and return that.
  tid_t insert_tid_string(const std::string& tid_name) {
    tid_t next_tid = static_cast<tid_t>(name_id_map.size());
    auto lookup = name_id_map.find(tid_name);
    if ( lookup == name_id_map.end()) {
      name_id_map.insert({tid_name, next_tid});
      tid_names.push_back(tid_name);
      return next_tid;
    } else {
      return lookup->second;
    }
  }

  // get the string name given an id
  const std::string& get_tid_name(tid_t id) {
    return (id < tid_names.size()) ? tid_names[id] : invalid_name;
  }


   // Add guide exon with TID and transcript start
  void addGuideExon(uint32_t start, uint32_t end, const tid_t& tid, //const std::string &tid,
                     char strand) {
     IntervalTree *tree = getTreeForStrand(strand);
     if (tree)
       tree->insert(start, end, tid);
   }
 
   // Call this after loading all transcript data
   // TODO: only build when first see TID in read
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
  IntervalNode *getNextNode(IntervalNode *node, const tid_t& tid, //) {const std::string &tid,
                             char strand) {
     IntervalTree *tree = getTreeForStrand(strand);
     return tree->getNextNodeForTid(node, tid);
   }
 
   // Get previous node in chain for a specific TID
  IntervalNode *getPrevNode(IntervalNode *node, const tid_t& tid, // const std::string &tid,
                             char strand) {
     IntervalTree *tree = getTreeForStrand(strand);
     return tree->getPrevNodeForTid(node, tid);
   }
};

struct MateInfo {
   uint32_t mate_index;
   tid_t transcript_id;
   uint32_t match_pos;
   uint32_t mate_size;
   bool same_transcript; // true if both mates map to same transcript
   bool valid_pair;      // true if this mate pair should be output
 
   MateInfo() : same_transcript(false), valid_pair(false) {}
};

struct ReadInfo {
  std::set<std::tuple<tid_t, uint32_t>> matches;
   GSamRecord *brec;
   bool valid_read;
   uint32_t read_index;
   uint32_t read_size;
 
   // Mate information
   std::string mate_name;
   bool is_paired;
   bool is_first_mate;
   bool is_reverse;
   bool mate_is_reverse;
   std::map<tid_t, MateInfo *> mate_info;
 
   ReadInfo()
       : matches(), brec(nullptr), valid_read(false), read_index(0),
         read_size(0), mate_name(), is_paired(false), is_first_mate(false),
         is_reverse(false), mate_is_reverse(false), mate_info() {}
 
   ~ReadInfo() {
     // brec is deleted by CReadAln
     for (auto &pair : mate_info) {
       auto info = std::get<1>(pair);
       delete info;
    }
  }
};

// Collect all refguide transcripts for a single genomic sequence
struct GRefData {
	GList<GffObj> rnas; 			// all transcripts on this genomic seq
	int gseq_id;
	const char* gseq_name;
	GRefData(int gid =- 1): rnas(false,true,false), gseq_id(gid), gseq_name(NULL) {
	gseq_id = gid;
	if (gseq_id >= 0)
		gseq_name = GffObj::names->gseqs.getName(gseq_id);
	}

	void add(GffReader* gffr, GffObj* t) {
		if (gseq_id >= 0) {
		if (gseq_id != t->gseq_id)
			GError("Error: invalid call to GRefData::add() - different genomic sequence!\n");
		}
		else { //adding first transcript, initialize storage
			gseq_id = t->gseq_id;
			gseq_name = t->getGSeqName();
			if (gffr->gseqtable[gseq_id] == NULL)
				GError("Error: invalid genomic sequence data (%s)!\n",gseq_name);
			rnas.setCapacity(gffr->gseqtable[gseq_id]->fcount);
		}
		rnas.Add(t);
		t->isUsed(true);
		//setLocus(t); //use the GRefLocus::mexons to quickly find an overlap with existing loci, or create a new one
	}

	bool operator==(GRefData& d) {
		return gseq_id == d.gseq_id;
	}
	bool operator<(GRefData& d) {
		return (gseq_id<d.gseq_id);
	}
};

// Bundle status
enum BundleStatus {
	BUNDLE_STATUS_CLEAR=0, 			// Available for loading/prepping
	BUNDLE_STATUS_LOADING=1, 		// Being prepared by the main thread (there can be only one)
	BUNDLE_STATUS_READY=2 			// Ready to be processed, or being processed
};

struct CBundle {
	int len;
	float cov;
	float multi;
	int start_node;  // id of start node in bundle of same strand
	int lastnodeid; // id of last node added to bundle
	CBundle(int _len = 0, float _cov = 0, float _multi = 0, int _start =- 1, int _last =- 1):
			len(_len), cov(_cov), multi(_multi), start_node(_start), lastnodeid(_last) {}
};

struct CGraphinfo {
	int ngraph;
	int nodeno;

	CGraphinfo(int ng = -1, int nnode = -1): ngraph(ng), nodeno(nnode) {}
};

struct CGJunc {
	int leftnode;
	int rightnode;
	double cov; // ngood
	double goodcov; // ngood_reads
	CGJunc(int n1=0,int n2=0,double _cov=0,double _goodcov=0):leftnode(n1),rightnode(n2),cov(_cov),goodcov(_goodcov){}
};


struct CGNode {
	int id;    // initial id in graphno
	bool last; // if this is last node (to be linked to sink later)
	bool keep; // if I keep it in the final count (true by default)
	bool merge; // if this node needs to be merged to its adjacent node
	bool future;
	CGNode(int _id=0,bool _last=false,bool _keep=true, bool _merge=false, bool _future=false):id(_id),last(_last),keep(_keep),merge(_merge),future(_future){}
};

struct GEdge { // guide edge
	// if val < endval then this is start; otherwise it is end
	uint val;  // value of the boundary
	uint endval; // value of the other exon boundary shared with val
	int strand;
	bool operator<(const GEdge& o) const {
		return(val<o.val || (val==o.val && strand<o.strand));
	}
	bool operator==(const GEdge& o) const {
		return(val==o.val && strand==o.strand);
	}
	GEdge(uint _val=0,uint _endval=0,int _strand=0):val(_val),endval(_endval),strand(_strand) {}
};

struct CGraphnode:public GSeg {
	int nodeid;
	float cov;
	float capacity; // sum of all transcripts abundances exiting and through node
	float rate; // conversion rate between in and out transfrags of node
	//float frag; // number of fragments included in node
	GVec<int> child;
	GVec<int> parent;
	GBitVec childpat;
	GBitVec parentpat;
	GVec<int> trf; // transfrags that pass the node
	bool hardstart:1; // verified/strong start
	bool hardend:1;	// verified/strong end
	//CGraphnode(int s=0,int e=0,unsigned int id=MAX_NODE,float nodecov=0,float cap=0,float r=0,float f=0):GSeg(s,e),nodeid(id),cov(nodecov),capacity(cap),rate(r),frag(f),child(),parent(),childpat(),parentpat(),trf(){}
	CGraphnode(int s=0,int e=0,unsigned int id=MAX_NODE,float nodecov=0,float cap=0,float r=0):GSeg(s,e),
			nodeid(id),cov(nodecov),capacity(cap),rate(r),child(),parent(),childpat(),parentpat(),trf(),hardstart(false),hardend(false){}
};

// # 0: strand; 1: start; 2: end; 3: nreads; 4: nreads_good;
struct CJunction:public GSeg {
	char strand; //-1,0,1
	char guide_match; //exact match of a ref intron?
	char consleft; // -1,0,1 -1 is not set up, 0 is non consensus, 1 is consensus
	char consright; // -1,0,1 -1 is not set up, 0 is non consensus, 1 is consensus
	double nreads;
	double nreads_good;
	double leftsupport;
	double rightsupport;
	double nm; // number of reads with a high nm (high mismatch)
	double mm; // number of reads that support a junction with both anchors bigger than longintronanchor
	CJunction(int s=0,int e=0, char _strand=0):GSeg(s,e),
			strand(_strand), guide_match(0), consleft(-1), consright(-1),nreads(0),nreads_good(0),
			leftsupport(0),rightsupport(0),nm(0),mm(0) {}
	bool operator<(CJunction& b) {
		if (start<b.start) return true;
		if (start>b.start) return false;
		if (end<b.end) return true;
		if (end>b.end) return false;
		if (strand>b.strand) return true;
		return false;
	}
	bool operator==(CJunction& b) {
		return (start==b.start && end==b.end && strand==b.strand);
	}
};


// Read Alignment
 struct CReadAln : public GSeg {
   // DEBUG ONLY:
   // GStr name;
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
   GPVec<CJunction> juncs;
   GSamRecord *brec; // store BAM record
 
   CReadAln(char _strand = 0, short int _nh = 0, int rstart = 0, int rend = 0)
       : GSeg(rstart, rend), // name(rname),
         strand(_strand), nh(_nh), len(0), read_count(0), unitig(false),
         longread(false), pair_count(), pair_idx(), mate_genomic_pos(0),
         mate_transcript_id(), mate_found(false), segs(), juncs(false), brec() {}
 
   CReadAln(CReadAln &rd) : GSeg(rd.start, rd.end) { // copy contructor
     strand = rd.strand;
     nh = rd.nh;
     len = rd.len;
     read_count = rd.read_count;
     unitig = rd.unitig;
     longread = rd.longread;
     pair_count = rd.pair_count;
     pair_idx = rd.pair_idx;
     juncs = rd.juncs;
   }
 
   int overlapSegLen(CReadAln *r) {
 
     if (r->start > end || start > r->end)
       return 0;
 
     int i = 0;
     int j = 0;
     int len = 0;
     while (i < segs.Count()) {
       if (segs[i].end < r->segs[j].start)
         i++;
       else if (r->segs[j].end < segs[i].start)
         j++;
       else { // there is overlap
         len += segs[i].overlapLen(r->segs[j].start, r->segs[j].end);
         if (segs[i].end < r->segs[j].end)
           i++;
         else
           j++;
       }
       if (j == r->segs.Count())
         break;
     }
     return len;
   }
   ~CReadAln() { delete brec; }

  friend std::ostream& operator<<(std::ostream&, const CReadAln&);
};

inline std::ostream& operator<<(std::ostream& os, const CReadAln& aln) {
  os << "[CReadAln](strand: " << aln.strand 
     << ", nh: " << aln.nh << ", len: " << aln.len
     << ", read_count: " << aln.read_count << ", unitig: " << aln.unitig
     << ", long_read: " << aln.longread 
     << ", mate_genomic_pos: " << aln.mate_genomic_pos 
     << ", mate_transcript_id: " << aln.mate_transcript_id
     << ", mate_found: " << aln.mate_found << ")";
  return os;
}

struct GReadAlnData {
	GSamRecord* brec;
	char strand; // -1, 0, 1
	int nh;
	int hi;
	GPVec<CJunction> juncs;

	GReadAlnData(GSamRecord* bamrec = nullptr, char nstrand = 0, int num_hits = 0,
	             int hit_idx = 0)
	    : brec(bamrec), strand(nstrand), nh(num_hits), hi(hit_idx), juncs(true) { }

	~GReadAlnData() = default;
};

struct RC_Feature { //exon or intron of a reference transcript

	uint id; 				//feature id (>0), +1 to the index either in bundle_RC_exons/introns
	GVec<uint> t_ids; 		//transcripts owning this feature
	 				  		// index in the BundleData::keepguides array + 1
	int l; int r; //genomic coordinates for the feature
	char strand;
	mutable uint rcount; //# read alignments overlapping this feature (>5bp overlaps for exons;
	                     // exact coord. match for introns)
	mutable uint ucount; //# uniquely mapped reads overlapping/matching this ref feature
	mutable double mrcount; //multi-map weighted read counts overlapping/matching this feature

	mutable double movlcount; //exons only: multi-map weighted sum of overlap lengths

    double avg;
    double stdev;
    double mavg;
    double mstdev;

	struct PCompare {
	 	bool operator()(const RC_Feature* p1, const RC_Feature* p2) {
	 		return (*p1 < *p2);
	 	}
	};

	RC_Feature(int l0 = 0, int r0 = 0, char s = '.', uint fid = 0, uint tid = 0): id(fid), t_ids(1), l(l0), r(r0),
		strand(s), rcount(0),ucount(0),mrcount(0), movlcount(0), avg(0), stdev(0), mavg(0), mstdev(0) {

		if (l > r) { int t = l; l = r; r = t; }
		if (tid > 0) t_ids.Add(tid);
	}

	RC_Feature(const RC_Feature& seg): id(seg.id), t_ids(seg.t_ids), l(seg.l), r(seg.r),
			strand(seg.strand), rcount(0), ucount(0), mrcount(0), movlcount(0), avg(0), stdev(0), mavg(0), mstdev(0) {
	}

	RC_Feature(const RC_Feature& seg, uint tid): id(seg.id), t_ids(1), l(seg.l), r(seg.r),
		strand(seg.strand), rcount(0), ucount(0), mrcount(0), movlcount(0), avg(0), stdev(0), mavg(0), mstdev(0) {
		if (l > r) { int t = l; l = r; r = t; }
		if (tid > 0) t_ids.Add(tid);
	}

	bool operator<(const RC_Feature& o) const {
		if (l != o.l) return (l < o.l);
     	if (r != o.r) return (r < o.r);
    	if (strand == '.' || o.strand == '.') return false;
    	if (strand != o.strand) return (strand < o.strand);
     return false;
	 }
	bool operator==(const RC_Feature& o) const {
	 //if (id == o.id) return true;
	 return (l==o.l && r==o.r &&
		 (strand == o.strand || strand == '.' || o.strand == '.'));
	 }
	bool strandCompatible(const RC_Feature& o) const {
		 return (strand == '.' || o.strand == '.' || strand == o.strand);
	}
	//WARNING: the overlap checks IGNORE strand!
	bool overlap(int hl, int hr) const {
	  if (hl>hr) { int t=hl; hl=hr; hr=t; }
      return (l<=hr && r<=hl);
	  }
	bool overlap(int hl, int hr, int minovl) const {
	  if (hl>hr) { int t=hl; hl=hr; hr=t; }
      hl+=minovl;hr-=minovl;
      return (l<=hr && r<=hl);
	  }
	uint ovlen(int hl, int hr) const {
     if (hl>hr) { int t=hl; hl=hr; hr=t; }
     if (l<hl) {
        if (hl>r) return 0;
        return (hr>r) ? r-hl+1 : hr-hl+1;
        }
       else { //hl<=l
        if (l>hr) return 0;
        return (hr<r)? hr-l+1 : r-l+1;
        }
	 }
};

struct RC_ExonOvl {
	RC_Feature* feature; //pointer to an item of RC_BundleData::g_exons
	int mate_ovl; // = 1 if the mate of this read overlaps the same exon
	int ovlen;
	bool operator<(const RC_ExonOvl& o) const {
		if (mate_ovl!=o.mate_ovl)
			return (mate_ovl>o.mate_ovl);
		if (ovlen!=o.ovlen)
			return (ovlen>o.ovlen);
		if (feature->r-feature->l != o.feature->r-o.feature->l)
			return (feature->r-feature->l > o.feature->r-o.feature->l);
		if (feature->strand != o.feature->strand)
			return (feature->strand<o.feature->strand);
		return (feature->l<o.feature->l);
	} //operator <
	bool operator==(const RC_ExonOvl& o) const {
		return (mate_ovl==o.mate_ovl && ovlen==o.ovlen && feature==o.feature);
	}
	RC_ExonOvl(RC_Feature* f=NULL, int olen=0, int movl=0):feature(f),
			mate_ovl(movl), ovlen(olen) {
	}
};

struct RC_TData { //storing RC data for a transcript
	//only used with -B (full Ballgown data)
	GffObj* ref_t;
	uint t_id;
	int l;
	int r;
	char in_bundle; // 1 if used by read bundles (present in keepguides),
	                // 2 if all introns are covered by at least one read, 3 if it is stored to be printed
	//GRefLocus* locus; //pointer to a locus info
	int eff_len;
	double cov;
	double fpkm;
	//char strand;
    GPVec<RC_Feature> t_exons;
    GPVec<RC_Feature> t_introns;

	void rc_addFeatures(uint& c_e_id, GList<RC_Feature>& fexons, GPVec<RC_Feature>& edata,
	                      uint& c_i_id, GList<RC_Feature>& fintrons, GPVec<RC_Feature>& idata) {
		GASSERT(ref_t);
  		GffObj& m = *(ref_t);
  		int ecache_idx = fexons.Count()-1;
  		int icache_idx = fintrons.Count()-1;
  		
  		for (int i = 0; i < m.exons.Count(); ++i)  {
    		addFeature((int)m.exons[i]->start, (int)m.exons[i]->end, t_exons, c_e_id, fexons, edata, ecache_idx);
    		//if (i==0) e_idx_cache=ecache_idx;
    		if (i>0) { //store intron
      		//if (i==1) i_idx_cache=icache_idx;
      		addFeature(m.exons[i-1]->end+1, m.exons[i]->start-1, t_introns, c_i_id,
    		  fintrons, idata, icache_idx);
    	} //for each intron
  	} //for each exon
						  }
	void addFeature(int fl, int fr, GPVec<RC_Feature>& fvec, uint& f_id,
			          GList<RC_Feature>& fset, GPVec<RC_Feature>& fdata,
					  int& cache_idx) {
		//f_id is the largest f_id inserted so far in fset
		bool add_new = true;
		RC_Feature* newseg=new RC_Feature(fl, fr, ref_t->strand, 0, this->t_id);
		//RC_Feature* newfeature=NULL;
		int fit=cache_idx<0 ? fset.Count()-1 : cache_idx;
		int fp_id=-1;
		if (fset.Count()>0) {
			if (*newseg < *(fset[fit])) {
			bool eq=false;
			while (*newseg < *(fset[fit]) || (eq = (*newseg==*(fset[fit])))) {
				if (eq) {
				add_new = false;
				fp_id = fset[fit]->id; //fset[fit]->id;
				break;
				}
				//newseg< fset[fit]
				--fit;
				if (fit<0) break; //newseg should be inserted at 0
			} //while newseg<fset[fit]
			if (add_new) ++fit;
				// newseg < fset[fit+1]
				//we'll insert newseg at position fit+1
			}
			else { //newseg >= *fset[fit]
			bool eq=false;
			while (*(fset[fit]) < *newseg || (eq = (*newseg==*(fset[fit])))) {
				if (eq) {
				add_new = false;
				fp_id = fset[fit]->id;
				break;
				}
				++fit;
				if (fit==fset.Count()) {
					//newseg should be appended to the list
					break;
				}
			}
			}
		} //check existing set
		if (add_new) { //did not see this feature before
			newseg->id = ++f_id;
			if (fit<0) fit = fset.Add(newseg);
			else fset.sortInsert(fit, newseg);
			if (fit<0) {
			GError("Error: feature %d-%d (%c) already in feature set!\n",
						newseg->l, newseg->r, newseg->strand);
			}
			cache_idx=fit;
			fp_id=fdata.Add(newseg)+1;

			GASSERT((uint)fdata.Count()==f_id);
		}
		else { //feature seen before, update its parent list
		fdata[fp_id-1]->t_ids.Add(this->t_id);
		delete newseg;
		}

		GASSERT(fdata[fp_id-1]->id==(uint)fp_id);
		fvec.Add(fdata[fp_id-1]);
	}

	RC_TData(GffObj& s, uint id=0):ref_t(&s), t_id(id), l(s.start), r(s.end),
			in_bundle(0), eff_len(s.covlen), cov(0), fpkm(0), //strand(s.strand),
			t_exons(false), t_introns(false) { //, e_idx_cache(-1), i_idx_cache(-1) {
	}

    bool operator<(const RC_TData& o) const {
    	if (l != o.l) return (l < o.l);
    	if (r != o.r) return (r < o.r);
    	if (char c=(ref_t->strand - o.ref_t->strand)) return (c<0);
	    return (strcmp(ref_t->getID(), o.ref_t->getID())<0);
    }
    bool operator==(const RC_TData& o) const {
    	if (t_id!=0 && o.t_id!=0) return (t_id==o.t_id);
    	return (l==o.l && r==o.r && ref_t->strand == o.ref_t->strand &&
    			strcmp(ref_t->getID(),o.ref_t->getID())==0);
    }
};

struct RC_BundleData {
	int init_lmin;
	int lmin;
	int rmax;
	GPVec<RC_TData> g_tdata; //raw counting data for all transcripts in this bundle
	// RC_TData::t_id-1 = the element index in this array
	GList<RC_Feature> g_exons; //set of guide exons in this bundle, sorted by start coordinate
	GList<RC_Feature> g_introns; //set of guide introns in this bundle, sorted by start coordinate
	//RC_FeatIt xcache; //cache the first exon overlapping xcache_pos to speed up exon-overlap queries (findOvlExons())
	int xcache; //exons index of the first exon overlapping xcache_pos
	int xcache_pos; // left coordinate of last cached exon overlap query (findOvlExons())

	// local (bundle) stable order tables of guide features
	GPVec<RC_TData>* bundle_RC_tdata; 				//pointer to the global RC tdata table
	// RC_Feature::id-1 = the index in these arrays:
	GPVec<RC_Feature>* bundle_RC_exons;  			// pointer to local exon RC data
	GPVec<RC_Feature>* bundle_RC_introns;			// pointer to local intron RC data
	//local exon/intron ids within the bundle
	uint c_exon_id;
	uint c_intron_id;
	//--
	RC_BundleData(int t_l=0, int t_r=0, GPVec<RC_TData>* rc_tdata=NULL,
		GPVec<RC_Feature>* rc_exons=NULL,GPVec<RC_Feature>* rc_introns=NULL):
		init_lmin(0), lmin(t_l), rmax(t_r),  g_tdata(false),
		// features:(sorted, free, unique)
		g_exons(true, false, true), g_introns(true, false, true),
		xcache(0), xcache_pos(0),
		bundle_RC_tdata(rc_tdata), bundle_RC_exons(rc_exons), bundle_RC_introns(rc_introns),
		c_exon_id(0), c_intron_id(0)
	{
		
		bundle_RC_exons  = new GPVec<RC_Feature>(true);
		bundle_RC_introns= new GPVec<RC_Feature>(true);
	}

	~RC_BundleData() {
		
		delete bundle_RC_exons;
		delete bundle_RC_introns;
	}

	uint addTranscript(GffObj& t) { //should return the guide index in *guides_RC_tdata
	
		if (lmin==0 || lmin>(int)t.start) lmin=t.start; 
		if (rmax==0 || rmax<(int)t.end) rmax=t.end; 

		GASSERT(t.uptr); //we should always have a RC_TData for each guide
		RC_TData* tdata=(RC_TData*)(t.uptr);
		g_tdata.Add(tdata);
		tdata->rc_addFeatures(c_exon_id, g_exons, *bundle_RC_exons,
			c_intron_id, g_introns, *bundle_RC_introns);
		return tdata->t_id;
	}

	bool findOvlExons(GArray<RC_ExonOvl>& exovls, int hl, int hr, char strand='.',
											int mate_pos=0, bool update_cache=true) {
		//exovls should be clear, unless the caller knows what s/he's doing
		bool hasOverlaps=false;
		if (g_exons.Count()==0) return false;
		RC_Feature q(hl, hr);
		int xstart=0;
		bool no_cache=(xcache_pos==0 || xcache_pos>hl);
		if (no_cache) {
			if (update_cache) {
				//xcache=exons.end();
				xcache=g_exons.Count()-1;
				xcache_pos=0;
			}
		}
		else xstart=xcache; //must have a valid value
		bool upd_cache(update_cache);
		int last_checked_exon=g_exons.Count()-1;
		for (int p=xstart;p < g_exons.Count();++p) {
			last_checked_exon=p;
			int l=g_exons[p]->l;
			int r=g_exons[p]->r;
			if (l > hr) break;
			if (hl > r) continue;
			//exon overlap here
			int ovlen=0;
			if (hl<l) {
				ovlen = ( hr<r ? hr-l+1 : r-l+1 );
			}
			else { // l<=hl
				ovlen= ( hr<r ? hr-hl+1 : r-hl+1 );
			}
			if (upd_cache) {
				//cache first overlap
				xcache=p;
				upd_cache=false;
			}
			if (strand!='.' && strand!=g_exons[p]->strand) continue; //non-matching strand
			int mate_ovl=0;
			if (mate_pos && mate_pos+10>l && mate_pos+5<r)
				mate_ovl=1; //mate read likely overlaps this exon
			if (mate_ovl || ovlen>=5) {
				//TODO: check this, arbitrary ovl minimum of 5bp
				hasOverlaps=true;
				RC_ExonOvl fovl(g_exons[p], ovlen, mate_ovl);
				exovls.Add(fovl);
			}
		}
		if (update_cache) {
			if (upd_cache) xcache=last_checked_exon; //there was no overlap found
			xcache_pos=hl;
		}
		return hasOverlaps;
	}
	
	RC_Feature* findIntron(int hl, int hr, char strand) {
	int fidx=0;
	RC_Feature* r=NULL;
	RC_Feature t(hl, hr, strand);
	if (g_introns.Found(&t, fidx))
		r=g_introns[fidx];
	return r;
	}
};

void rc_update_exons(RC_BundleData& rc);


// Bundle data structure, holds input data parsed from BAM file
struct BundleData {
	BundleStatus status;
	
	int idx; 							// index in the main bundles array
	int start;
	int end;
	unsigned long num_reads; 			// number of reads in this bundle
	double num_fragments; 				// aligned read/pairs

	GStr refseq; 						// reference sequence name
	char* gseq; 						// actual genomic sequence for the bundle
	GList<CReadAln> reads;			    // list of reads in bundle
	GList<CJunction> junction;
	GPVec<GffObj> guides;
	RC_BundleData* rc_data;

	BundleData():status(BundleStatus::BUNDLE_STATUS_CLEAR), idx(0), start(0), end(0),
			num_reads(0), num_fragments(0), refseq(), gseq(NULL), reads(false, true), // readlist(sorted, free_elements)
			junction(true, true, true), guides(false), rc_data(NULL) { }

	// This is only called when the bundle is valid and ready to be processed
	void getReady(int curr_start, int curr_end) {
		start = curr_start;
		end = curr_end;

		// Tag all these guides
		for (int i = 0; i < this->guides.Count(); ++i) {
			RC_TData* tdata = (RC_TData*)(guides[i]->uptr);
			tdata->in_bundle = 1;
		}
		status = BundleStatus::BUNDLE_STATUS_READY;
	}

	// Create new RC_BundleData
	void rc_init(GffObj* t, GPVec<RC_TData>* rc_tdata,
			     GPVec<RC_Feature>* rc_edata, GPVec<RC_Feature>* rc_idata) {

		if (rc_data == NULL) {
			rc_data = new RC_BundleData(t->start, t->end, rc_tdata, rc_edata, rc_idata);
		}
	}
	
	void keepGuide(GffObj* t, GPVec<RC_TData>* rc_tdata=NULL,
			GPVec<RC_Feature>* rc_edata=NULL, GPVec<RC_Feature>* rc_idata=NULL) {
		if (rc_data==NULL) {
			rc_init(t, rc_tdata, rc_edata, rc_idata);
		}
		guides.Add(t);
		t->udata=(int)rc_data->addTranscript(*t);
	}

	bool evalReadAln(GReadAlnData& alndata, char& strand);

	void Clear() {

		guides.Clear();
		reads.Clear();
		junction.Clear();
		
		start = 0;
		end = 0;
		status = BundleStatus::BUNDLE_STATUS_CLEAR;
		num_reads = 0;
		num_fragments = 0;

		delete rc_data;
		GFREE(gseq);
		rc_data = NULL;
	}

	~BundleData() {
		Clear();
	}
};

extern GStr tmp_path;
extern GStr cram_ref;
extern bool keepTempFiles;

struct BamIO {
	protected:
		GSamReader* reader = nullptr;
		GSamWriter* writer = nullptr;
		GSamRecord* current_record = nullptr;
		GStr input_bam;
		GStr output_bam;
		GStr header_sam;

	public:
		BamIO(const GStr& in_bam, const GStr& out_bam, const GStr& sam_header)
			: input_bam(in_bam), output_bam(out_bam), header_sam(sam_header) {}

		void start() {
			reader = new GSamReader(input_bam.chars());
			writer = new GSamWriter(output_bam.chars(), header_sam.chars());

			GSamRecord* next_rec = reader->next();
             if (next_rec) {
                current_record = next_rec;
            }
		}

		GSamRecord* next() {
			GSamRecord* result = current_record;

            current_record = reader->next();  // null at EOF
			
            return result;
		}

		void write(GSamRecord* rec) {
			if (writer != nullptr && rec != nullptr) {
				writer->write(rec);
			}
		}

		int32_t get_tid(const char* transcript_id) {
    		return writer->get_tid(transcript_id);
		}

		void stop() {
			
			if (reader != nullptr) {
				delete reader;
				reader = nullptr;
			}
			if (writer != nullptr) {
				delete writer;
				writer = nullptr;
			}
			if (current_record != nullptr) {
				delete current_record;
				current_record = nullptr;
			}
		}
};

struct WorkerArgs {
    GPVec<BundleData>* bundle_queue;
    BamIO* io;

	WorkerArgs(GPVec<BundleData>* q, BamIO* d) : bundle_queue(q), io(d) {}
};

#endif
