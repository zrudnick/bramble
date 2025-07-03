#ifndef MAIN_H
#define MAIN_H

#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <memory>
#include <unordered_map>

#include "GArgs.h"
#include "GStr.h"
#include "gff.h"
#include "GSam.h"
#include "GBitVec.h"
#include "time.h"
#include "GHashMap.hh"

using namespace std;
 
#ifndef NOTHREADS				// Use threads
#define THREADING_ENABLED 1
#else
#define THREADING_ENABLED 0
#endif

#define MAX_NODE 1000000
#define SMALL_EXON 35 			// Exons smaller than this have a tendency to be missed by long read data
#define RUNOFF_DIST 200		    // Reads at what distance should be considered part of separate bundles?
#define LONG_INTRON_ANCHOR 25
#define MISMATCH_FRAC 0.02

struct IntervalNode {
    uint start, end;
    uint max_end;  // Maximum end value in subtree (interval tree augmentation)
	uint height;
    std::set<std::string> tids;
    std::map<std::string, uint> tid_cum_len;
    IntervalNode* left;
    IntervalNode* right;

	// Maps from TID to previous/next nodes in genomic order for that TID
    std::map<std::string, IntervalNode*> tid_prev_nodes;  // TID -> previous node with same TID
    std::map<std::string, IntervalNode*> tid_next_nodes;  // TID -> next node with same TID
    
    IntervalNode(uint s, uint e) : start(s), end(e), max_end(e), height(1),
				left(nullptr), right(nullptr) {}
};

class IntervalTree {
private:
    IntervalNode* root;
	std::set<std::string> all_tids;

	int getHeight(IntervalNode* node) {
        return node ? node->height : 0;
    }
    
    int getBalance(IntervalNode* node) {
        return node ? getHeight(node->left) - getHeight(node->right) : 0;
    }
    
    void updateHeight(IntervalNode* node) {
        if (node) {
			uint left_height = getHeight(node->left);
			uint right_height = getHeight(node->right);
            node->height = 1 + std::max(left_height, right_height);
        }
    }
    
    void updateMaxEnd(IntervalNode* node) {
        if (!node) return;
        
        node->max_end = node->end;
        if (node->left) {
            node->max_end = std::max(node->max_end, node->left->max_end);
        }
        if (node->right) {
            node->max_end = std::max(node->max_end, node->right->max_end);
        }
    }

	IntervalNode* rotateRight(IntervalNode* y) {
        IntervalNode* x = y->left;
        IntervalNode* T2 = x->right;
        
        // Rotation
        x->right = y;
        y->left = T2;
        
        updateHeight(y);
        updateHeight(x);
        updateMaxEnd(y);
        updateMaxEnd(x);
        
        return x;
    }
    
    IntervalNode* rotateLeft(IntervalNode* x) {
        IntervalNode* y = x->right;
        IntervalNode* T2 = y->left;
        
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
    
    IntervalNode* getParent(IntervalNode* node) {
        if (!node || node == root) return nullptr;
        return findParent(root, node);
    }
    
    IntervalNode* findParent(IntervalNode* current, IntervalNode* target) {
        if (!current || current == target) return nullptr;
        if (current->left == target || current->right == target) return current;
        
        IntervalNode* leftResult = findParent(current->left, target);
        if (leftResult) return leftResult;
        return findParent(current->right, target);
    }
    
    void destroyTree(IntervalNode* node) {
        if (node) {
            destroyTree(node->left);
            destroyTree(node->right);
            delete node;
        }
    }
    
    void inorderTraversal(IntervalNode* node, std::vector<IntervalNode*>& result) {
        if (node) {
            inorderTraversal(node->left, result);
            result.push_back(node);
            inorderTraversal(node->right, result);
        }
    }

	IntervalNode* insertNodeBalanced(IntervalNode* node, IntervalNode* newNode) {
        if (!node) return newNode;
        
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
            (newNode->start == node->left->start && newNode->end < node->left->end))) {
            return rotateRight(node);
        }
        
        // Case 2: Right Right
        if (balance < -1 && (newNode->start > node->right->start || 
            (newNode->start == node->right->start && newNode->end > node->right->end))) {
            return rotateLeft(node);
        }
        
        // Case 3: Left Right
        if (balance > 1 && (newNode->start > node->left->start || 
            (newNode->start == node->left->start && newNode->end > node->left->end))) {
            node->left = rotateLeft(node->left);
            return rotateRight(node);
        }
        
        // Case 4: Right Left
        if (balance < -1 && (newNode->start < node->right->start || 
            (newNode->start == node->right->start && newNode->end < node->right->end))) {
            node->right = rotateRight(node->right);
            return rotateLeft(node);
        }
        
        return node;
    }
    
public:
    IntervalTree() : root(nullptr) {}

    ~IntervalTree() {
        destroyTree(root);
    }

	void insert(uint start, uint end, const std::string& tid) {
        // Insert tid to set of all TIDs
        all_tids.insert(tid);

        // Find tree nodes that overlap
        auto overlapping = findAllOverlapping(start, end);

        // No overlaps: create new node for entire range
        if (overlapping.empty()) {
            auto newNode = new IntervalNode(start, end);
            newNode->tids.insert(tid);
            //newNode->height = 1;
            insertNode(newNode);
            return;
        }

        // Check for exact match
        for (auto& node : overlapping) {
            if (node->start == start && node->end == end) {
                node->tids.insert(tid);
                return;
            }
        }

        // *^*^*^ *^*^*^ *^*^*^ *^*^*^
    	// Split nodes
    	// *^*^*^ *^*^*^ *^*^*^ *^*^*^

        std::set<uint> split_points;
        for (auto& node : overlapping) {
            if (node->start < start && start < node->end) {
                split_points.insert(start);
            }
            if (node->start < end && end < node->end) {
                split_points.insert(end);
            }
        }

        // Perform all splits in one pass
        for (uint split_point : split_points) {
            auto nodes_at_point = findAllOverlapping(split_point, split_point + 1);
            for (auto& node : nodes_at_point) {
                if (node->start < split_point && node->end > split_point) {
                    splitInterval(node, split_point);
                }
            }
        }

        // Now find all intervals that should contain this TID
        auto final_overlapping = findAllOverlapping(start, end);

        // Track which parts of [start, end) are covered by existing nodes
        std::vector<std::pair<uint, uint>> covered_ranges;
        covered_ranges.reserve(final_overlapping.size());

        for (auto& node : final_overlapping) {
            uint overlap_start = std::max(node->start, start);
            uint overlap_end = std::min(node->end, end);

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
        std::vector<IntervalNode*> new_nodes;
        
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
        for (auto* node : new_nodes) {
            node->tids.insert(tid);
            //node->height = 1;
            insertNode(node);
        }
    }
    
private:

	void insertNode(IntervalNode* newNode) {
        root = insertNodeBalanced(root, newNode);
    }
    
    std::vector<IntervalNode*> findAllOverlapping(uint start, uint end) {
        std::vector<IntervalNode*> result;
        findOverlappingHelper(root, start, end, result);
        return result;
    }
    
    void findOverlappingHelper(IntervalNode* node, uint start, uint end, 
                              std::vector<IntervalNode*>& result) {
        if (!node) return;
        
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
	std::vector<std::pair<uint, uint>> mergeRanges(std::vector<std::pair<uint, uint>>& ranges) {
		if (ranges.empty()) return ranges;
		
		std::vector<std::pair<uint, uint>> merged;
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
	void splitInterval(IntervalNode* node, uint split_point) {
		std::set<std::string> tids = node->tids;
		
		// Create new node for second half
		auto newNode = new IntervalNode(split_point, node->end);
		newNode->tids = tids;
		// newNode->height = 1;
		
		node->end = split_point; 		// update original node to represent first half
		insertNode(newNode); 			// insert the new node into tree
	}

	void getTidNodesHelper(IntervalNode* node, const std::string& tid, 
							std::vector<IntervalNode*>& result) {
		if (!node) return;
		
		if (node->tids.count(tid)) {
			result.push_back(node);
		}
		
		getTidNodesHelper(node->left, tid, result);
		getTidNodesHelper(node->right, tid, result);
	}

	// Get all nodes containing a specific TID
	std::vector<IntervalNode*> getTidNodes(const std::string& tid) {
		std::vector<IntervalNode*> result;
		getTidNodesHelper(root, tid, result);
		return result;
	}

	// Build chain for a specific TID
	void buildTidChain(const std::string& tid) {
		// Get all nodes with this TID in genomic order
		auto tidNodes = getTidNodes(tid);
		
		// Sort by genomic position
		std::sort(tidNodes.begin(), tidNodes.end(),
				[](IntervalNode* a, IntervalNode* b) {
					if (a->start != b->start) return a->start < b->start;
					return a->end < b->end;
				});
		
		// Build bidirectional chain
		for (size_t i = 0; i < tidNodes.size(); ++i) {
			if (i > 0) {
				tidNodes[i]->tid_prev_nodes[tid] = tidNodes[i-1];
				tidNodes[i-1]->tid_next_nodes[tid] = tidNodes[i];
			}
		}
	}
    
public:
	// For debugging (print_tree)
    std::vector<IntervalNode*> getOrderedIntervals() {
        std::vector<IntervalNode*> result;
        inorderTraversal(root, result);
        return result;
    }
    
    // Find all intervals that overlap with the given range
    std::vector<IntervalNode*> findOverlapping(uint start, uint end) {
        std::vector<IntervalNode*> result;
        findOverlappingHelper(root, start, end, result);
        std::sort(result.begin(), result.end(),
                 [](IntervalNode* a, IntervalNode* b) 
				 { return a->start < b->start; });
        return result;
    }
    
    // Get next node in chain for a specific TID
    IntervalNode* getNextNodeForTid(IntervalNode* node, const std::string& tid) {
        auto it = node->tid_next_nodes.find(tid);
        return (it != node->tid_next_nodes.end()) ? it->second : nullptr;
    }
    
    // Get previous node in chain for a specific TID
    IntervalNode* getPrevNodeForTid(IntervalNode* node, const std::string& tid) {
        auto it = node->tid_prev_nodes.find(tid);
        return (it != node->tid_prev_nodes.end()) ? it->second : nullptr;
    }

	// Build all TID chains after tree construction is complete
	void buildAllTidChains() {
		for (const auto& tid : all_tids) {
			buildTidChain(tid);
		}
	}

	// Precompute cumulative lengths so that match positions can be calculated
    void precomputeCumulativeLengths(const std::string& tid) {
        auto tidChain = getTidNodes(tid);

		// Need to sort by genomic position
		// TODO: store these chains so they don't need to be resorted
		std::sort(tidChain.begin(), tidChain.end(),
				[](IntervalNode* a, IntervalNode* b) {
					if (a->start != b->start) return a->start < b->start;
					return a->end < b->end;
				});
        
        uint cumulative = 0;
		uint prev_end = -1;
        for (size_t i = 0; i < tidChain.size(); ++i) {
            auto* node = tidChain[i];
            
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
        for (const auto& tid : all_tids) {
            precomputeCumulativeLengths(tid);
        }
    }

	// Check for cumulative length
	uint findCumulativeLength(IntervalNode* node, const std::string& tid) {
        auto it = node->tid_cum_len.find(tid);
        return (it != node->tid_cum_len.end()) ? it->second : 0;
    }
    
};

// g2t tree using interval tree
struct g2tTree {
    IntervalTree* fw_tree;		// tree for forward strand guides
	IntervalTree* rc_tree;      // tree for reverse strand guides
    
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
    
    // Add guide exon with TID and transcript start
    void addGuideExon(uint start, uint end, const std::string& tid, char strand) {
		IntervalTree* tree = getTreeForStrand(strand);
		if (tree) tree->insert(start, end, tid);
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
    std::vector<IntervalNode*> getIntervals(uint readStart, uint readEnd, char strand) {
		IntervalTree* tree = getTreeForStrand(strand);
        if (!tree) return std::vector<IntervalNode*>();

        return tree->findOverlapping(readStart, readEnd);
    }

	// Get cumulative previous size of exons from transcript
	uint getCumulativeLength(IntervalNode* node, const std::string& tid, char strand) {
        IntervalTree* tree = getTreeForStrand(strand);
		return tree->findCumulativeLength(node, tid);
    }

	// Get next node in chain for a specific TID
    IntervalNode* getNextNode(IntervalNode* node, const std::string& tid, char strand) {
        IntervalTree* tree = getTreeForStrand(strand);
        return tree->getNextNodeForTid(node, tid);
    }
    
    // Get previous node in chain for a specific TID
    IntervalNode* getPrevNode(IntervalNode* node, const std::string& tid, char strand) {
        IntervalTree* tree = getTreeForStrand(strand);
        return tree->getPrevNodeForTid(node, tid);
    }

};

struct MateInfo {
    uint mate_index;
	std::string transcript_id;
	uint match_pos;
	uint mate_size;
	bool same_transcript; 		// true if both mates map to same transcript
	bool valid_pair;			// true if this mate pair should be output
    
    MateInfo() : same_transcript(false), valid_pair(false) {}
};

struct ReadInfo {
	std::set<std::tuple<std::string, uint>> matches;
	GSamRecord* brec;
	bool valid_read;
	uint read_index;
	uint read_size;
	uint nh_i;
	bool wrote_alr;		// this read (a mate) was already written to bam
	
	// Mate information
	std::string mate_name;
	bool is_paired;
	bool is_first_mate;
	bool is_reverse;
	bool mate_is_reverse;
	std::map<std::string, MateInfo*> mate_info;

	ReadInfo() : matches(), brec(nullptr), valid_read(false), read_index(0), read_size(0), nh_i(0),
	             wrote_alr(false), mate_name(), is_paired(false), is_first_mate(false), 
				 is_reverse(false), mate_is_reverse(false), mate_info() {}
	
    ~ReadInfo() {
		// brec is deleted by CReadAln
		for (auto& pair : mate_info) {
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
	BUNDLE_STATUS_CLEAR = 0, 			// Available for loading/prepping
	BUNDLE_STATUS_LOADING = 1, 			// Being prepared by the main thread (there can be only one)
	BUNDLE_STATUS_READY = 2 			// Ready to be processed, or being processed
};

// Read Alignment
struct CReadAln:public GSeg {
	char strand;                // 1, 0 (unknown), -1 (reverse)
	short int nh;
	uint len;
	float read_count;           // keeps count for all reads (including paired and unpaired)
	bool unitig:1;			    // set if read come from an unitig
	bool longread:1;	        // set if read comes from long read data
	GVec<float> pair_count;     // keeps count for all paired reads
	GVec<int> pair_idx;         // keeps indeces for all pairs in assembly mode, or all reads that were collapsed in merge mode
	
	// For mate position calculations
	uint mate_genomic_pos;
	std::string mate_transcript_id;
	bool mate_found;

	GVec<GSeg> segs;            //"exons"					
	GSamRecord* brec;			// store BAM record

	CReadAln(char _strand = 0, short int _nh = 0,
			int rstart = 0, int rend = 0)
			: GSeg(rstart, rend), //name(rname),
			strand(_strand), nh(_nh), len(0), 
			read_count(0), unitig(false), 
			longread(false), pair_count(), 
			pair_idx(), mate_genomic_pos(0),
			mate_transcript_id(), mate_found(false),
			segs(), brec() { }
			
	CReadAln(CReadAln &rd):GSeg(rd.start,rd.end) { // copy contructor
		strand = rd.strand;
		nh = rd.nh;
		len = rd.len;
		read_count = rd.read_count;
		unitig = rd.unitig;
		longread = rd.longread;
		pair_count = rd.pair_count;
		pair_idx = rd.pair_idx;
	}

	~CReadAln() { 
		delete brec;
	}
};


struct GroupKey {
    uint start;
    std::string cigar_str;
    
    // for std::map
    bool operator<(const GroupKey& other) const {
        if (start != other.start) return start < other.start;
        return cigar_str < other.cigar_str;
    }
    
    bool operator==(const GroupKey& other) const {
        return start == other.start && cigar_str == other.cigar_str;
    }
};

struct AlnGroups {
    std::map<GroupKey, uint> key_to_group;
    std::vector<std::vector<uint32_t>> groups;
    
    AlnGroups() : key_to_group(), groups() {}

    std::string getCigar(bam1_t* b) {
        uint32_t n_cigar = b->core.n_cigar;
        uint32_t* cigar = bam_get_cigar(b);
        
        std::string cigar_str;
        cigar_str.reserve(n_cigar * 8);
        
        for (uint32_t i = 0; i < n_cigar; i++) {
            uint32_t op_len = bam_cigar_oplen(cigar[i]);
            char op_char = bam_cigar_opchr(cigar[i]);
            cigar_str += std::to_string(op_len) + op_char;
        }
        return cigar_str;
    }

    void Add(CReadAln* read, uint n) {
        if (!read || !read->brec) return;
        bam1_t* b = read->brec->get_b();
        
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
    
    ~AlnGroups() {
        Clear();
    }
};

// Bundle data structure, holds input data parsed from BAM file
struct BundleData {
	BundleStatus status;				// bundle status
	int idx; 							// index in the main bundles array
	int start;							// bundle start
	int end;							// bundle end

	GStr refseq; 						// reference sequence name
	char* gseq; 						// genomic sequence for the bundle
	GList<CReadAln> reads;			    // all reads in bundle
	GPVec<GffObj> guides;				// all guides in bundle

	BundleData():status(BundleStatus::BUNDLE_STATUS_CLEAR), idx(0), start(0), end(0), 
			refseq(), gseq(NULL), reads(false, true), guides(false) //rc_data(NULL) { } // sorted, free elements
			{}

	// Called when the bundle is valid and ready to be processed
	void getReady(int curr_start, int curr_end) {
		start = curr_start;
		end = curr_end;
		status = BundleStatus::BUNDLE_STATUS_READY;
	}

	void keepGuide(GffObj* t) {
		guides.Add(t);
	}

	void Clear() {
		guides.Clear();
		reads.Clear();
		GFREE(gseq);

		start = 0;
		end = 0;
		status = BundleStatus::BUNDLE_STATUS_CLEAR;
	}

	~BundleData() {
		Clear();
	}
};

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