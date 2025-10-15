
#include <algorithm>
#include <memory>
#include <set>
#include <string>
#include <vector>
#include <unordered_set>

#include "types.h"
#include "bundles.h"
#include "bramble.h"
#include "g2t.h"
#include "evaluate.h"
#include "bam.h"

#ifndef NOTHREADS
#include "GThreads.h"
#endif

extern GFastMutex bam_io_mutex;  // protects BAM io

extern bool LONG_READS;
extern bool USE_FASTA;
extern bool SOFT_CLIPS;
extern bool STRICT;
const uint32_t overhang_threshold = 8;

namespace bramble {

  ReadEvaluator::~ReadEvaluator() {}

  // Evaluate read group and return matches
  ReadEvaluationResult 
  ReadEvaluator::evaluate(BundleData *bundle, uint read_index, 
                          g2tTree *g2t) {
    return ReadEvaluationResult{};
  }

  std::string ReadEvaluator::reverse_complement(const std::string &seq) {
    std::string rc;
    rc.reserve(seq.length());

    for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
      switch (*it) {
      case 'A': rc += 'T'; break;
      case 'T': rc += 'A'; break;
      case 'G': rc += 'C'; break;
      case 'C': rc += 'G'; break;
      case 'N': rc += 'N'; break;
      default: rc += *it; break;
      }
    }
    return rc;
  }

  // @requires: gseq != null
  std::string ReadEvaluator::extract_sequence(char *gseq, uint start, 
                                              uint length, char strand) {
    std::string seq(gseq + start, length);
    if (strand == '-') seq = reverse_complement(seq);
    return seq;
  }

  std::tuple<uint32_t, uint32_t> 
  ReadEvaluator::get_exon_coordinates(BundleData *bundle, GSeg exon, char strand) {
    uint32_t exon_start, exon_end;
    if (strand == '+') {
      exon_start = exon.start - bundle->start;
      exon_end = exon.end - bundle->start;
    } else {
      exon_start = bundle->end - exon.end;
      exon_end = bundle->end - exon.start;
    }
    return std::make_tuple(exon_start, exon_end);
  }
    
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
  pos_t ReadEvaluator::get_match_pos(IntervalNode *interval, 
                                    tid_t tid, char read_strand,
                                    g2tTree *g2t, uint exon_start,
                                    bool used_backwards_overhang) {
    IntervalNode *first_interval = interval;
     
    if (used_backwards_overhang) {
      first_interval = g2t->getPrevNode(interval, tid, read_strand);
    }

    uint32_t interval_start = first_interval->start;

    uint32_t prev_node_sum =
        g2t->getCumulativeLength(first_interval, tid, read_strand);
    uint32_t match_start = (exon_start > interval_start) ? 
      (exon_start - interval_start) : 0;
    return (pos_t)(match_start + prev_node_sum);
  }

  /**
   * Check if we can match forwards to another interval 
   * (USE_FASTA mode)
   *
   */
  bool ReadEvaluator::search_forward(IntervalNode *interval, uint exon_end,
                                    char read_strand,
                                    std::set<tid_t> &exon_tids, g2tTree *g2t,
                                    BundleData *bundle) {

    std::set<tid_t> tmp_exon_tids;
    uint32_t overhang_length = exon_end - interval->end;
    if (overhang_length > overhang_threshold) return false;

    uint32_t overhang_start;
    if (read_strand == '+') {
      overhang_start = interval->end;
    } else {
      uint32_t relative_end = bundle->end - bundle->start;
      overhang_start = relative_end - interval->end - overhang_length;
    }

    // Extract overhang sequence from genomic sequence
    if (bundle->gseq == nullptr) return false;
    std::string overhang_seq = extract_sequence(bundle->gseq, overhang_start,
      overhang_length, read_strand);

    for (const auto &tid : exon_tids) {
      IntervalNode *next_interval = g2t->getNextNode(interval, tid, read_strand);

      if (next_interval && next_interval->start >= overhang_start) {
        uint32_t check_length = 
          std::min(overhang_length, next_interval->end - next_interval->start);
        
        uint32_t next_overhang_start;
        if (read_strand == '+' || read_strand == 1) {
          next_overhang_start = next_interval->start;
        } else {
          next_overhang_start =
            (bundle->end - next_interval->start - check_length - bundle->start);
        }

        std::string guide_seq = extract_sequence(bundle->gseq, 
          next_overhang_start, check_length, read_strand);

        if (overhang_seq == guide_seq) {
          tmp_exon_tids.insert(tid);
        }
      }
    }
    
    std::swap(exon_tids, tmp_exon_tids);
    tmp_exon_tids.clear();

    return !exon_tids.empty();
  }

  /**
   * Check if we can match backwards to another interval 
   * (USE_FASTA mode)
   *
   */
  bool ReadEvaluator::search_backward(IntervalNode *interval, uint exon_start,
                                      char read_strand,
                                      std::set<tid_t> &exon_tids, g2tTree *g2t,
                                      BundleData *bundle) {

    std::set<tid_t> tmp_exon_tids;
    uint32_t overhang_length = interval->start - exon_start;
    if (overhang_length > overhang_threshold) return false;

    uint32_t overhang_start;
    if (read_strand == '+') {
      overhang_start = exon_start;
    } else {
      uint32_t relative_end = bundle->end - bundle->start;
      overhang_start = relative_end - exon_start - overhang_length;
    }

    // Extract overhang sequence from genomic sequence
    if (bundle->gseq == nullptr) return false;
    std::string overhang_seq = extract_sequence(bundle->gseq, overhang_start,
      overhang_length, read_strand);

    for (const auto &tid : exon_tids) {
      IntervalNode *prev_interval = g2t->getPrevNode(interval, tid, read_strand);

      if (prev_interval &&
          prev_interval->end <= overhang_start + overhang_length) {
        uint32_t check_length =
            std::min(overhang_length, prev_interval->end - prev_interval->start);
        
        uint32_t prev_overhang_start;
        if (read_strand == '+') {
          prev_overhang_start = prev_interval->end - check_length;
        } else {
          prev_overhang_start = bundle->end - prev_interval->end - bundle->start;
        }

        std::string guide_seq = extract_sequence(
          bundle->gseq, prev_overhang_start, check_length, read_strand);

        if (overhang_seq == guide_seq) {
          tmp_exon_tids.insert(tid);
        }
      }
    }
    std::swap(exon_tids, tmp_exon_tids);
    tmp_exon_tids.clear();

    return !exon_tids.empty();
  }

  std::vector<char> 
  ReadEvaluator::get_strands_to_check(CReadAln* read) {
    if (read->strand == '+') return {'+'}; 
    if (read->strand == '-') return {'-'}; 
    return {'+', '-'};  // unstranded: try forward first, then reverse
  }

  // Dynamic boundary tolerance
  uint32_t ReadEvaluator::boundary_tolerance(uint32_t exon_len) {
    if (STRICT) return 0;
    // 5% of exon length, max 500bp
    if (LONG_READS) return 500; //std::min(exon_len/2, 500u);
    // 50% of exon length, max 20bp
    else return std::min(exon_len/2, 20u);
  }

  std::set<tid_t> 
  ReadEvaluator::get_candidate_tids(std::vector<IntervalNode *> intervals) {
    if (intervals.empty()) return {};
    
    std::set<tid_t> tids(intervals[0]->tids.begin(), intervals[0]->tids.end());
    
    for (size_t i = 1; i < intervals.size(); ++i) {
      
      // Check for continuity in intervals
      if (intervals[i]->start != intervals[i-1]->end) {
        return {};  // gap found
      }

      // todo: allow gaps
      // requires individual cigars
      
      // Keep only TIDs present in all intervals
      auto it = tids.begin();
      while (it != tids.end()) {
        if (std::find(intervals[i]->tids.begin(), intervals[i]->tids.end(), *it) 
            == intervals[i]->tids.end()) {
          it = tids.erase(it);
        } else {
          ++it;
        }
      }
    }
    
    return tids;
  }

  void ReadEvaluator::hide_small_exon(Cigar& cigar,
                                      uint32_t exon_start, uint32_t exon_end,
                                      bool is_first_exon, bool is_last_exon) {
    uint32_t match_size = exon_end - exon_start + 1;

    // First exon: soft clip whole thing
    if (is_first_exon && SOFT_CLIPS) {
      cigar.add_operation(match_size, BAM_CSOFT_CLIP);
    }
    
    // Middle exon: entire thing is an insertion
    if (!is_first_exon && !is_last_exon) {
      cigar.add_operation(match_size, BAM_CINS);
    }
    
    // Last exon: soft clip whole thing
    if (is_last_exon && SOFT_CLIPS) {
      cigar.add_operation(match_size, BAM_CSOFT_CLIP);
    }
  }

  bool ReadEvaluator::check_first_exon(uint32_t exon_start, uint32_t interval_start,
                                      uint32_t exon_end, uint32_t interval_end,
                                      bool is_last_exon, uint32_t max_clip_size) {
    
    uint32_t tolerance = boundary_tolerance(exon_end - exon_start);

    if (!SOFT_CLIPS) {
      if (is_last_exon && 
        (exon_start < interval_start || exon_end >= interval_end)) {
        return false;
      } else if (exon_start < interval_start || exon_end - tolerance >= interval_end) {
        return false;
      }
    } else {
      if (exon_end - tolerance >= interval_end) {
        return false;
      } else if (exon_start < interval_start) {
        uint32_t clip_size = interval_start - exon_start;
        if ((clip_size >= max_clip_size)){
          return false;
        }
      }
    }
    return true;
  }
  
  void ReadEvaluator::filter_tids(IntervalNode* first_interval,
                                  IntervalNode* prev_last_interval,
                                  std::set<tid_t> &candidate_tids,
                                  uint32_t exon_start,
                                  g2tTree* g2t, char strand) {
    // Keep only TIDs with correct continuity
    auto it = candidate_tids.begin();
    while (it != candidate_tids.end()) {
      auto tid = *it;
      
      // Check if this TID appears in the first interval
      if (std::find(first_interval->tids.begin(), 
                    first_interval->tids.end(), 
                    tid) == first_interval->tids.end()) {
        it = candidate_tids.erase(it);
        continue;
      }
      
      // Check continuity: previous exon for this TID must match prev_last_interval
      // todo: allow insertions
      // requires individual cigars
      auto tid_last_interval = g2t->getPrevNode(first_interval, tid, strand);
      if (tid_last_interval == prev_last_interval) {
        ++it;
      } else {
        it = candidate_tids.erase(it);
      }
    }
  }

  bool ReadEvaluator::check_middle_exon(uint32_t exon_start, 
                                        uint32_t interval_start,
                                        uint32_t exon_end, 
                                        uint32_t interval_end,
                                        IntervalNode* first_interval,
                                        IntervalNode* prev_last_interval,
                                        std::set<tid_t> &candidate_tids,
                                        g2tTree* g2t, char strand) {

    uint32_t tolerance = boundary_tolerance(exon_end - exon_start);

    if ((exon_start + tolerance < interval_start) || 
      (exon_end - tolerance >= interval_end)) {
      return false;
    }

    filter_tids(first_interval, prev_last_interval, 
      candidate_tids, exon_start, g2t, strand);
    return true;
  }

  bool ReadEvaluator::check_last_exon(uint32_t exon_start, 
                                      uint32_t interval_start,
                                      uint32_t exon_end, 
                                      uint32_t interval_end,
                                      IntervalNode* first_interval,
                                      IntervalNode* prev_last_interval,
                                      std::set<tid_t> &candidate_tids,
                                      g2tTree* g2t, char strand, 
                                      uint32_t max_clip_size) {

    uint32_t tolerance = boundary_tolerance(exon_end - exon_start);

    if (!SOFT_CLIPS) {
      if ((exon_start + tolerance < interval_start) || 
        (exon_end >= interval_end)) {
        return false;
      }
    } else {
      if (exon_start + tolerance < interval_start) {
        return false;
      } else if (exon_end >= interval_end) {
        uint32_t clip_size = exon_end - (interval_end - 1);
        if ((clip_size >= max_clip_size)){
          return false;
        }
      }
    }

    filter_tids(first_interval, prev_last_interval, 
      candidate_tids, exon_start, g2t, strand);
    return true;
  }

  // Calculate similarity between exon chain and transcript intervals
  double 
  ReadEvaluator::exon_chain_similarity(GVec<GSeg>& read_exons,
                                      const std::vector<std::vector<IntervalNode*>>& 
                                      all_intervals,
                                      tid_t tid, char strand, BundleData* bundle) {
    uint32_t total_overlap = 0;
    uint32_t total_read_length = 0;
    
    for (int i = 0; i < read_exons.Count(); i++) {
      GSeg exon = read_exons[i];
      auto [exon_start, exon_end] = get_exon_coordinates(bundle, exon, strand);

      total_read_length += (exon_end - exon_start);
      
      // Check intervals that were found for this exon
      for (const auto& interval : all_intervals[i]) {
        if (std::find(interval->tids.begin(), interval->tids.end(), tid) != 
            interval->tids.end()) {
          uint32_t overlap_start = std::max(exon_start, interval->start);
          uint32_t overlap_end = std::min(exon_end, interval->end);
          if (overlap_end > overlap_start) {
            total_overlap += (overlap_end - overlap_start);
          }
        }
      }
    }
    
    return (total_read_length > 0) ? ((double)total_overlap / total_read_length) : 0.0;
  }

  void ReadEvaluator::build_exon_cigar(Cigar& cigar, bool is_first_exon, 
                                      bool is_last_exon, 
                                      uint32_t exon_start, uint32_t exon_end,
                                      uint32_t interval_start, uint32_t interval_end) {

    // Handle start boundary
    if (exon_start < interval_start) {
      uint32_t overhang = interval_start - exon_start;
      if (is_first_exon && SOFT_CLIPS) cigar.add_operation(overhang, BAM_CSOFT_CLIP);
      else if (!is_first_exon) cigar.add_operation(overhang, BAM_CINS);
    }
    
    // Add match for overlapping region
    uint32_t overlap_start = std::max(exon_start, interval_start);
                                           // interval_end - 1 to make it inclusive
    uint32_t overlap_end = std::min(exon_end, interval_end - 1);  
    if (overlap_end >= overlap_start) {
      uint32_t match_length = overlap_end - overlap_start + 1;
      cigar.add_operation(match_length, BAM_CMATCH);
    }
    
    // Handle end boundary
    if (exon_end >= interval_end) {
      uint32_t overhang = exon_end - interval_end + 1;
      if (is_last_exon && SOFT_CLIPS) cigar.add_operation(overhang, BAM_CSOFT_CLIP);
      else if (!is_last_exon) cigar.add_operation(overhang, BAM_CINS);
    }
  }

  bool ReadEvaluator::compile_matches(std::vector<ExonChainMatch> &matches,
                                      std::set<tid_t> candidate_tids, 
                                      IntervalNode *first_interval, 
                                      char strand, g2tTree* g2t,
                                      uint32_t exon_start, 
                                      bool is_first_exon,
                                      bool used_backwards_overhang) {
    if (is_first_exon) {
      // Initialize matches with TIDs and positions
      for (const auto &tid : candidate_tids) {
        uint32_t match_pos = get_match_pos(first_interval, 
          tid, strand, g2t, exon_start, used_backwards_overhang);
        ExonChainMatch match;
        match.tid = tid;
        match.pos = match_pos;
        match.gaps_count = 0;
        match.total_gap_size = 0;
        match.total_coverage = 0;
        match.soft_clip_front = 0;
        match.soft_clip_back = 0;
        matches.emplace_back(match);
      }

    } else {
      auto it = matches.begin();
      while (it != matches.end()) {
        if (candidate_tids.count(it->tid)) {
          ++it;
        } else {
          it = matches.erase(it);
        }
      }
      if (matches.empty()) {
        return false;
      }
    }
  return true;
  }

  void ReadEvaluator::filter_by_similarity(std::vector<ExonChainMatch> &matches,
                                           std::vector<std::vector<IntervalNode*>> 
                                           all_matched_intervals,
                                           GVec<GSeg> read_exons, char strand,
                                           BundleData *bundle, 
                                           double similarity_threshold) {
    auto it = matches.begin();
    while (it != matches.end()) {
      double similarity = exon_chain_similarity(read_exons, 
        all_matched_intervals, it->tid, strand, bundle);
      if (similarity > similarity_threshold) {
        ++it;
      } else {
        it = matches.erase(it);
      }
    }
  }
 
  ReadEvaluationResult 
  ReadEvaluator::evaluate_exon_chains(BundleData *bundle, read_id_t id, 
                                      g2tTree *g2t, 
                                      ReadEvaluationConfig config) {
    CReadAln *read = bundle->reads[id];
    GVec<GSeg> read_exons = read->segs;
    uint32_t exon_count = read_exons.Count();

    IntervalNode *first_interval = nullptr;
    IntervalNode *last_interval = nullptr;
    IntervalNode *prev_last_interval = nullptr;
    uint32_t interval_start;
    uint32_t interval_end;
    
    std::vector<ExonChainMatch> matches;
    std::vector<std::vector<IntervalNode*>> all_matched_intervals;
    Cigar cigar;  // single CIGAR for all matches

    auto strands_to_check = get_strands_to_check(read);
    for (char strand : strands_to_check) {
      matches.clear(); 
      cigar.clear();
      all_matched_intervals.clear();
      all_matched_intervals.resize(exon_count);
      bool strand_failed = false;
      uint32_t total_exon_matches = exon_count;

      for (uint j = 0; j < exon_count; j++) {
        auto [exon_start, exon_end] = get_exon_coordinates(bundle, 
          read_exons[j], strand);
        bool is_first_exon = (j == 0);
        bool is_last_exon = (j == exon_count - 1);

        // Find all guide intervals that contain the read exon
        auto intervals = g2t->getIntervals(exon_start, exon_end, strand);
        if (intervals.empty()) {
          if (config.ignore_small_exons 
            && (exon_end - exon_start <= config.small_exon_size)) {
            hide_small_exon(cigar, exon_start, exon_end, 
              is_first_exon, is_last_exon); 
            total_exon_matches--;
            continue;
          }
          strand_failed = true;
          break;
        }

        // Get candidate tids for this exon
        auto candidate_tids = get_candidate_tids(intervals);

        // Store intervals for this exon
        all_matched_intervals[j] = intervals;
        prev_last_interval = last_interval;
        first_interval = intervals[0];
        last_interval = intervals[intervals.size() - 1];
        interval_start = first_interval->start;
        interval_end = last_interval->end;

        bool used_backwards_overhang = false;

        // Boundary checks depending on current exon
        if (is_first_exon) {
          strand_failed = !(check_first_exon(exon_start, interval_start,
            exon_end, interval_end, is_last_exon, config.max_clip_size));
        } else if (is_last_exon) {
           strand_failed = !(check_last_exon(exon_start, interval_start, 
            exon_end, interval_end, first_interval, prev_last_interval, 
            candidate_tids, g2t, strand, config.max_clip_size));
        } else {
          strand_failed = !(check_middle_exon(exon_start, interval_start, 
            exon_end, interval_end, first_interval, prev_last_interval, 
            candidate_tids, g2t, strand));
        }
        if (strand_failed) break;

        // Build CIGAR as we process each exon
        build_exon_cigar(cigar, is_first_exon, is_last_exon, exon_start, 
          exon_end, interval_start, interval_end);

        // Update matches
        strand_failed = !compile_matches(matches, candidate_tids, 
          first_interval, strand, g2t, exon_start, is_first_exon, 
          used_backwards_overhang);
        if (strand_failed) break;

        // Fail if we run out of candidiate TIDs
        if (candidate_tids.empty()) {
          strand_failed = true;
          break;
        }

        // After the last exon, filter matches by similarity
        if (config.filter_similarity && is_last_exon) 
          filter_by_similarity(matches, all_matched_intervals, read_exons,
            strand, bundle, config.similarity_threshold);
      }

      // If this strand worked, return matches
      if (!strand_failed && !matches.empty()) {
        if (total_exon_matches > 0) return {true, matches, strand, cigar};
        else return config.default_result;
      }
      
      // If this was a stranded read and it failed, give up
      if (read->strand != '.' && matches.empty()) {
        return config.default_result;
      }

      // If this strand failed, try the other strand
    }

    // If this was an unstranded read and both strands fail, 
    // give up
    return config.default_result;
  }

  ReadEvaluationResult 
  ShortReadEvaluator::evaluate(BundleData *bundle, read_id_t id, g2tTree *g2t) {
    ReadEvaluationConfig config = {
      20,                             // max clip size
      false,                          // filter similarity
      0.90,                           // similarity threshold
      true,                           // ignore small exons
      10,                             // small exon size
      {false, {}, '.', {}}            // default result
    }; 

    ReadEvaluationResult res = evaluate_exon_chains(bundle, id, g2t, config);
    return res;
  }

  ReadEvaluationResult 
  LongReadEvaluator::evaluate(BundleData *bundle, read_id_t id, g2tTree *g2t) {
    ReadEvaluationConfig config = {
      500,                            // max clip size
      false,                          // filter similarity
      0.60,                           // similarity threshold
      true,                           // ignore_small_exons
      500,                            // small exon size
      {false, {}, '.', {}}            // default result
    };   

    ReadEvaluationResult res = evaluate_exon_chains(bundle, id, g2t, config); 
    return res;
  }

} // namespace bramble