
#include <algorithm>
#include <memory>
#include <set>
#include <string>
#include <vector>
#include <unordered_set>
#include <cmath>
#include <random>

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

extern uint32_t dropped_reads;
extern uint32_t unresolved_reads;

std::string read_name;
double similarity_threshold;  // constant, record for bam.cpp

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
      exon_end = (exon.end + 1) - bundle->start;
    } else {
      exon_start = bundle->end - (exon.end + 1);
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
                                    tid_t tid, char strand,
                                    g2tTree *g2t, uint exon_start,
                                    bool used_backwards_overhang) {
    IntervalNode *first_interval = interval;
     
    if (used_backwards_overhang) {
      first_interval = g2t->getPrevNode(interval, tid, strand);
    }

    uint32_t interval_start = first_interval->start;
    pos_t prev_node_sum =
      g2t->getCumulativeLength(first_interval, tid, strand);
    std::string tid_string = g2t->getTidName(tid);
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

  std::set<tid_t> 
  ReadEvaluator::get_candidate_tids(std::vector<IntervalNode *> intervals,
                                    std::vector<ExonChainMatch> &matches,
                                    std::unordered_map<tid_t, std::shared_ptr<Cigar>> &subcigars,
                                    uint32_t exon_start, uint32_t exon_end,
                                    bool is_first_exon, bool is_last_exon, g2tTree *g2t) {
    if (intervals.empty()) 
      return {};

    if (STRICT || !LONG_READS) {
      std::set<tid_t> tids(intervals[0]->tids.begin(), intervals[0]->tids.end());

      for (uint32_t i = 1; i < intervals.size(); ++i) {
      
        // Check for continuity in intervals
        if (intervals[i]->start != intervals[i-1]->end) {
          return {};  // gap found
        }

        // Keep only TIDs present in all intervals
        auto it = tids.begin();
        while (it != tids.end()) {
          if (!intervals[i]->tids.count(*it)) {
            it = tids.erase(it);
          } else {
            ++it;
          }
        }
      }
      return tids;
    
    // LONG READS
    } else {
       // Add ALL possible tids from every interval
      std::set<tid_t> tids;
      uint32_t max_ins = 10;
      for (uint32_t i = 0; i < intervals.size(); ++i) {

        for (auto &tid_here : intervals[i]->tids) {
          tids.emplace(tid_here);
        }

        uint32_t match_length = intervals[i]->end - intervals[i]->start;
        if (match_length > max_ins) {
          auto it = tids.begin();
          while (it != tids.end()) {
            // If long interval, remove any TIDs not in this interval
            if (!intervals[i]->tids.count(*it)) {
              it = tids.erase(it);
            } else {
              ++it;
            }
          }
        }
      }

      // Determine TID-based cigar portions for this exon
      for (auto &tid : tids) {
        subcigars[tid] = std::make_shared<Cigar>();
      }

      auto not_in = [&](uint32_t start, uint32_t end, tid_t tid) {
        if (start >= intervals.size()) return true;
        end = std::min<uint32_t>(end, intervals.size());
        while (start < end) {
          if (intervals[start]->tids.count(tid)) return false;
          start++;
        }
        return true;
      };

      uint32_t n_intervals = intervals.size();
      for (uint32_t i = 0; i < n_intervals; i++) {

        // If intervals are non-contiguous, this represents insertions
        // Small insertions are common in long-read data
        if (i != 0 && intervals[i]->start != intervals[i-1]->end) {
          uint32_t insertion_length = intervals[i]->start - intervals[i-1]->end;
          if (insertion_length > max_ins) return {};
          for (auto &tid_here : intervals[i]->tids) {
            if (tids.count(tid_here)) 
              subcigars[tid_here]->add_operation(insertion_length, BAM_CINS);
          }
        }

        uint32_t overlap_start = std::max(exon_start, intervals[i]->start);
        uint32_t overlap_end = std::min(exon_end, intervals[i]->end); 
        if (overlap_end <= overlap_start) continue;

        uint32_t match_length = overlap_end - overlap_start;
        for (auto &tid : tids) {
          // Match
          if (intervals[i]->tids.count(tid)) {
            subcigars[tid]->add_operation(match_length, BAM_CMATCH);

          // Insertion / Soft Clip
          } else {
            if (is_first_exon && not_in(0, i, tid)) {
              subcigars[tid]->add_operation(match_length, BAM_CSOFT_CLIP);
            } else if (is_last_exon && not_in(i+1, n_intervals, tid)) {
              subcigars[tid]->add_operation(match_length, BAM_CSOFT_CLIP);
            } else {
              subcigars[tid]->add_operation(match_length, BAM_CINS);
            }
          }
        }
      }
      return tids;
    }
  }

  // Add operation to every CIGAR in matches
  void ReadEvaluator::add_operation_for_all(std::vector<ExonChainMatch> &matches,
                                            uint32_t length, uint8_t op, g2tTree *g2t, 
                                            char strand) {
    auto it = matches.begin();
    while (it != matches.end()) {
      std::shared_ptr<Cigar> match_cigar = it->align.cigar;
      match_cigar->add_operation(length, op);
      if (op != BAM_CSOFT_CLIP && op != BAM_CINS) {
        it->ref_consumed += length;
        it->total_coverage += length;
      }
      ++it;
    }
  };

  void ReadEvaluator::hide_small_exon(std::vector<ExonChainMatch> &matches,
                                      uint32_t exon_start, uint32_t exon_end,
                                      bool is_first_exon, bool is_last_exon, g2tTree *g2t, char strand) {

    uint32_t match_size = exon_end - exon_start;

    // First exon: soft clip whole thing
    if (is_first_exon && SOFT_CLIPS) {
      add_operation_for_all(matches, match_size, BAM_CSOFT_CLIP, g2t, strand);
    }
    
    // Middle exon: entire thing is an insertion
    if (!is_first_exon && !is_last_exon) {
      add_operation_for_all(matches, match_size, BAM_CINS, g2t, strand);
    }
    
    // Last exon: soft clip whole thing
    if (is_last_exon && SOFT_CLIPS) {
      add_operation_for_all(matches, match_size, BAM_CSOFT_CLIP, g2t, strand);
    }
  }
  
  // @requires: not called on first exon
  void ReadEvaluator::filter_tids(std::vector<IntervalNode *> intervals,
                                  std::unordered_map<tid_t, IntervalNode *> 
                                  prev_intervals,
                                  std::set<tid_t> &candidate_tids,
                                  uint32_t exon_start,
                                  g2tTree* g2t, char strand,
                                  std::unordered_map<tid_t, uint32_t> &gaps) {
    // Keep only TIDs with correct continuity
    auto it_c = candidate_tids.begin();
    while (it_c != candidate_tids.end()) {
      auto tid = *it_c;
      
      // Get first interval containing current TID
      IntervalNode* tmp = nullptr;
      if (strand == '+') {
        for (auto &interval : intervals) {
          if (interval->tids.count(tid)) {
            tmp = interval;
            break;
          }
        }
      } else { // "-"
        for (auto it = intervals.rbegin(); it != intervals.rend(); ++it) {
          auto &interval = *it;
          if (interval->tids.count(tid)) {
            tmp = interval;
            break;
          }
        }
      }
      
      // Shouldn't happen
      if (tmp == nullptr) {
        it_c = candidate_tids.erase(it_c);
        continue;
      }
      
      // Check continuity: previous exon for this TID must match prev_last_interval
      // don't allow gaps, it tanks accuracy
      // or allow gaps of like <10
      IntervalNode* tid_last_interval;
      if (strand == '+') tid_last_interval = g2t->getPrevNode(tmp, tid, strand);
      else tid_last_interval = g2t->getNextNode(tmp, tid, strand);

      // Long reads often have short deletions ("gaps") in reference to the reference
      uint32_t gap_size = 0;
      uint32_t max_gap = 10;
      
      if (LONG_READS && tid_last_interval != nullptr) {
        while (prev_intervals.count(tid) && tid_last_interval != prev_intervals[tid]) {
          gap_size += (tid_last_interval->end - tid_last_interval->start);
          if (gap_size > max_gap) break; // give up if gap if too large

          if (strand == '+') tid_last_interval = g2t->getPrevNode(tmp, tid, strand);
          else tid_last_interval = g2t->getNextNode(tmp, tid, strand);
          if (tid_last_interval == nullptr) break; // give up if we reach last possible node
        }
      }
      if (prev_intervals.count(tid) && 
        tid_last_interval == prev_intervals[tid]) {
        if (LONG_READS) gaps[tid] = gap_size;
        ++it_c;
      } else {
        it_c = candidate_tids.erase(it_c);
      }
    }
  }

  bool ReadEvaluator::check_first_exon(uint32_t exon_start, uint32_t interval_start,
                                      uint32_t exon_end, uint32_t interval_end,
                                      bool is_last_exon, uint32_t max_clip_size,
                                      uint32_t tolerance) {
    
    if (STRICT) tolerance = 0;

    if (!SOFT_CLIPS) {
      if (is_last_exon && 
        (exon_start < interval_start || exon_end > interval_end)) {
        return false;
      } else if (exon_start < interval_start || exon_end - tolerance > interval_end) {
        return false;
      }
    } else {
      if (exon_end - tolerance > interval_end) {
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

  bool ReadEvaluator::check_middle_exon(uint32_t exon_start, 
                                        uint32_t interval_start,
                                        uint32_t exon_end, 
                                        uint32_t interval_end,
                                        uint32_t tolerance) {

    if (STRICT) tolerance = 0;

    if ((exon_start + tolerance < interval_start) || 
      (exon_end - tolerance > interval_end)) {
      return false;
    }
    return true;
  }

  bool ReadEvaluator::check_last_exon(uint32_t exon_start, 
                                      uint32_t interval_start,
                                      uint32_t exon_end, 
                                      uint32_t interval_end,
                                      uint32_t max_clip_size,
                                      uint32_t tolerance) {

    if (STRICT) tolerance = 0;

    if (!SOFT_CLIPS) {
      if ((exon_start + tolerance < interval_start) || 
        (exon_end > interval_end)) {
        return false;
      }
    } else {
      if (exon_start + tolerance < interval_start) {
        return false;
      } else if (exon_end > interval_end) {
        uint32_t clip_size = exon_end - interval_end;
        if ((clip_size >= max_clip_size)){
          return false;
        }
      }
    }
    return true;
  }

  void ReadEvaluator::build_exon_cigar(std::vector<ExonChainMatch> &matches,
                                      std::unordered_map<tid_t, std::shared_ptr<Cigar>> &subcigars,
                                      bool is_first_exon, bool is_last_exon, 
                                      uint32_t exon_start, uint32_t exon_end,
                                      uint32_t interval_start, uint32_t interval_end,
                                      std::unordered_map<tid_t, uint32_t> &gaps,
                                      g2tTree* g2t, char strand) {
    // Long reads: Handle gaps between previous and current exon
    if (LONG_READS) {
      for (auto &match : matches) {
        tid_t tid = match.tid;
        std::shared_ptr<Cigar> cigar = match.align.cigar;
        if (gaps.count(tid) && gaps[tid] > 0) cigar->add_operation(gaps[tid], BAM_CDEL);
      }
    }

    // Handle start boundary
    if (exon_start < interval_start) {
      uint32_t overhang = interval_start - exon_start;
      if (is_first_exon && SOFT_CLIPS) add_operation_for_all(matches, overhang, BAM_CSOFT_CLIP, g2t, strand);
      else if (!is_first_exon) add_operation_for_all(matches, overhang, BAM_CINS, g2t, strand);
    }
    
    // Add match for overlapping region
    if (STRICT || !LONG_READS) {
      uint32_t overlap_start = std::max(exon_start, interval_start);
      uint32_t overlap_end = std::min(exon_end, interval_end);  
      if (overlap_end >= overlap_start) {
        uint32_t match_length = overlap_end - overlap_start;
        add_operation_for_all(matches, match_length, BAM_CMATCH, g2t, strand);
      }
    } else {
      auto it_m = matches.begin();
      while (it_m != matches.end()) {
        tid_t tid = it_m->tid;
        std::shared_ptr<Cigar> cigar = it_m->align.cigar;
        // Might happen if TID in previous exons
        // was not found for this exon
        auto it_s = subcigars.find(tid);
        if (it_s == subcigars.end()) {
          ++it_m;
          continue;
        }
        // If tid in subcigar, add subcigar to match cigar
        std::shared_ptr<Cigar> subcigar = it_s->second;
        for (uint32_t i = 0; i < subcigar->cigar.size(); i++) {
          uint32_t len = subcigar->cigar[i].first;
          uint8_t op = subcigar->cigar[i].second;
          cigar->add_operation(len, op);

          // Update coverage for similarity calculation
          if (op == BAM_CMATCH || op == BAM_CDEL) {
            if (op == BAM_CMATCH) it_m->total_coverage += len;
            it_m->ref_consumed += len;
          }
        }
        
        ++it_m;
      }
    } 
    
    // Handle end boundary
    if (exon_end > interval_end) {
      uint32_t overhang = exon_end - interval_end;
      if (is_last_exon && SOFT_CLIPS) add_operation_for_all(matches, overhang, BAM_CSOFT_CLIP, g2t, strand);
      else if (!is_last_exon) add_operation_for_all(matches, overhang, BAM_CINS, g2t, strand);
    }
  }

  void ReadEvaluator::compile_matches(std::vector<ExonChainMatch> &matches,
                                      std::set<tid_t> candidate_tids, 
                                      std::unordered_map<tid_t, IntervalNode *> &first_intervals, 
                                      char strand, g2tTree* g2t,
                                      uint32_t exon_start, uint32_t exon_end,
                                      uint32_t cumulative_length, uint32_t coverage_augment,
                                      bool used_backwards_overhang) {
    // Initialize matches with TIDs and positions
    for (const auto &tid : candidate_tids) {
      bool match_exists = false;
      auto it = matches.begin();
      while (it != matches.end()) {
        if (it->tid == tid) {
          match_exists = true;
          break;
        }
        ++it;
      }

      if (!match_exists) {
        ExonChainMatch match;
        match.tid = tid;
        pos_t pos = get_match_pos(first_intervals[tid], 
          tid, strand, g2t, exon_start, used_backwards_overhang);
        match.align.pos = pos;
        match.align.strand = strand;
        match.align.cigar = std::make_shared<Cigar>();
        match.transcript_length = g2t->getTranscriptLength(tid, strand);
        match.ref_consumed = 0;
        match.total_coverage = 0 + coverage_augment; // same as ref consumed
        match.total_read_length = cumulative_length;
          // current exon length added in update_matches
        match.align.similarity_score = 0;
        matches.emplace_back(match);

        if (cumulative_length > 0) {
          match.align.cigar->add_operation(cumulative_length, BAM_CSOFT_CLIP);
        }
      }
    }
  }

  bool ReadEvaluator::update_matches(std::vector<ExonChainMatch> &matches,
                                      std::set<tid_t> candidate_tids, 
                                      std::unordered_map<tid_t, IntervalNode *> &first_intervals, 
                                      char strand, g2tTree* g2t,
                                      uint32_t exon_start, uint32_t exon_end,
                                      bool is_first_exon,
                                      std::unordered_map<tid_t, uint32_t> &gaps,
                                      bool used_backwards_overhang) {
    
    auto it = matches.begin();
    while (it != matches.end()) {
      it->total_read_length += (exon_end - exon_start);

      if (!is_first_exon) {
        // For reads on negative strand, need to update
        // match position at each new exon
        tid_t tid = it->tid;
        if (strand == '-' && first_intervals.count(tid)) {
          pos_t pos = get_match_pos(first_intervals[tid], 
            tid, strand, g2t, exon_start, used_backwards_overhang);
          it->align.pos = pos;
        }

        // Remove tids unsupported by all exons
        if (candidate_tids.count(it->tid)) {
          ++it;
        } else {
          it = matches.erase(it);
        }
      } else {
        ++it;
      }
      
    }
    return !matches.empty();
  }

  void ReadEvaluator::get_prev_intervals(std::vector<IntervalNode *> intervals,
                                        std::set<tid_t> candidate_tids,
                                        std::unordered_map<tid_t, IntervalNode *>
                                        &prev_intervals, char strand) {
    // Update last seen intervals for each candidate tids
    // Previous candidate tids not in this interval set 
    // will simply not be updated
    for (tid_t tid : candidate_tids) {
      if (strand == '+') {
        for (int32_t i = intervals.size() - 1; i >= 0; i--) {
          IntervalNode* interval = intervals[i];
          if (interval->tids.count(tid)) {
            prev_intervals[tid] = interval;
            break;
          }
        }
      } else { // "-"
        for (int32_t i = 0; i < intervals.size(); i++) {
          IntervalNode* interval = intervals[i];
          if (interval->tids.count(tid)) {
            prev_intervals[tid] = interval;
            break;
          }
        }
      }
      
    }
  }

  void ReadEvaluator::get_first_intervals(std::vector<IntervalNode *> intervals,
                                          std::set<tid_t> candidate_tids,
                                          char strand,
                                          std::unordered_map<tid_t, IntervalNode *>
                                          &first_intervals) {
    for (tid_t tid : candidate_tids) {
        for (uint32_t i = 0; i < intervals.size(); i++) {
          IntervalNode* interval = intervals[i];
          if (interval->tids.count(tid)) {
            first_intervals[tid] = interval;
            break;
          }
        }

      // same for '+' and '-'
    }
  }

  int32_t get_rand(uint32_t x) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, x - 1);
    return dis(gen);
  }

  void ReadEvaluator::filter_by_similarity(std::vector<ExonChainMatch> &matches, 
                                           double similarity_threshold) {

    auto it = matches.begin();
    while (it != matches.end()) {
      uint32_t total_coverage = it->total_coverage;
      uint32_t total_read_length = it->total_read_length;
      double similarity = (total_read_length > 0) ? 
        ((double)total_coverage / total_read_length) : 0.0;
      if (similarity > similarity_threshold) {
        it->align.similarity_score = similarity;
        ++it;
      } else {
        it = matches.erase(it);
      }
    }

    if (matches.empty()) return;

    ExonChainMatch *best_match;
    tid_t best_score = -1;
    int32_t hit_index = 1;
    bool has_best = false;
    for (auto &match : matches) {
      match.align.hit_index = hit_index;
      hit_index++;
      if (best_score == -1) {
          best_score = match.align.similarity_score;
          best_match = &match;
      } else if (best_score != -1 && match.align.similarity_score > best_score) {
        best_score = match.align.similarity_score;
        best_match = &match;
        has_best = true;
      }
    }
    if (has_best) {
      best_match->align.primary_alignment = true;
    } else {
      int32_t rand_idx = get_rand(matches.size());
      matches[rand_idx].align.primary_alignment = true;
    }
    
  }

  // calculate FP and TP for long reads
  // venn diagrams
 
  ReadEvaluationResult 
  ReadEvaluator::evaluate_exon_chains(BundleData *bundle, read_id_t id, 
                                      g2tTree *g2t, 
                                      ReadEvaluationConfig config) {
    CReadAln *read = bundle->reads[id];
    GVec<GSeg> read_exons = read->segs;
    uint32_t exon_count = read_exons.Count();

    IntervalNode *first_interval = nullptr;
    IntervalNode *last_interval = nullptr;
    uint32_t interval_start;
    uint32_t interval_end;

    std::vector<ExonChainMatch> matches;
    std::vector<ExonChainMatch> matches_by_strand;
    std::unordered_map<tid_t, IntervalNode *> prev_intervals;

    auto strands_to_check = get_strands_to_check(read);
    for (char strand : strands_to_check) {
      matches.clear();
      prev_intervals.clear();
      bool strand_failed = false;
      uint32_t total_exon_matches = exon_count;
      uint32_t cumulative_length = 0;
      uint32_t coverage_augment = 0;

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
            hide_small_exon(matches, exon_start, exon_end, 
              is_first_exon, is_last_exon, g2t, strand); 
            cumulative_length += (exon_end - exon_start);
            total_exon_matches--;
            coverage_augment += (exon_end - exon_start);
            continue;
          }
          dropped_reads++;
          strand_failed = true;
          break;
        }

        // Get candidate tids for this exon
        std::unordered_map<tid_t, std::shared_ptr<Cigar>> subcigars;
        auto candidate_tids = get_candidate_tids(intervals, matches, subcigars,
          exon_start, exon_end, is_first_exon, is_last_exon, g2t);
        if (candidate_tids.empty()) {
          strand_failed = true;
          unresolved_reads++;
          break;
        }

        // Store intervals for this exon
        first_interval = intervals[0];
        last_interval = intervals[intervals.size() - 1];
        interval_start = first_interval->start;
        interval_end = last_interval->end;

        bool used_backwards_overhang = false;
        std::unordered_map<tid_t, uint32_t> gaps;

        // Boundary checks depending on current exon
        if (is_first_exon) {
          strand_failed = !(check_first_exon(exon_start, interval_start,
            exon_end, interval_end, is_last_exon, config.max_clip_size, 
            config.boundary_tolerance));
        } else if (is_last_exon) {
           strand_failed = !(check_last_exon(exon_start, interval_start, 
            exon_end, interval_end, config.max_clip_size, config.boundary_tolerance));
        } else {
          strand_failed = !(check_middle_exon(exon_start, interval_start, 
            exon_end, interval_end, config.boundary_tolerance));
        }
        if (strand_failed) {
          unresolved_reads++;
          break;
        }

        // Filter tids
        if (!is_first_exon) filter_tids(intervals, prev_intervals, 
          candidate_tids, exon_start, g2t, strand, gaps);
        if (candidate_tids.empty()) {
          strand_failed = true;
          unresolved_reads++;
          break;
        }

        std::unordered_map<tid_t, IntervalNode *> first_intervals;
        get_first_intervals(intervals, candidate_tids, strand, first_intervals);

        // Create exon chain matches 
        if (is_first_exon) {
          compile_matches(matches, candidate_tids, first_intervals, 
            strand, g2t, exon_start, exon_end, cumulative_length, 
            coverage_augment, used_backwards_overhang);
        }

        // Build CIGAR as we process each exon
        build_exon_cigar(matches, subcigars,
          is_first_exon, is_last_exon, exon_start, 
          exon_end, interval_start, interval_end, gaps, g2t, strand);

        strand_failed = !update_matches(matches, 
          candidate_tids, first_intervals, strand, g2t, 
          exon_start, exon_end, is_first_exon, gaps, 
          used_backwards_overhang);
        if (strand_failed) {
          unresolved_reads++;
          break;
        }

        // Fail if we run out of candidate TIDs
        if (candidate_tids.empty()) {
          unresolved_reads++;
          strand_failed = true;
          break;
        }

        if (!is_last_exon) get_prev_intervals(intervals, 
          candidate_tids, prev_intervals, strand);

        // After the last exon, filter matches by similarity
        if (config.filter_similarity && is_last_exon) 
          filter_by_similarity(matches, config.similarity_threshold);

        if (is_last_exon && matches.empty()) {
          strand_failed = true;
          break;
        }

        // Calculate at end for new matches in next exon
        cumulative_length += (exon_end - exon_start);
        if (total_exon_matches == 0) {
          unresolved_reads++;
          strand_failed = true;
          break;
        }

      } // end of exon loop

      // Handle short reads
      if (!LONG_READS) {
        // Worked on this strand, so we are done
        if (!strand_failed && !matches.empty()) {
          return {true, matches};
        }
        // Stranded read failed -> fail
        if (read->strand != '.') {
          return config.default_result;
        }
        // Unstranded: try the opposite strand next
        continue;
      }

      // Handle long reads
      if (LONG_READS) {
        // Stranded long reads: either success or failure
        if (read->strand != '.') {
          if (!strand_failed && !matches.empty()) {
            return {true, matches};
          } else {
            return config.default_result;
          }
        }

        // Unstranded long reads: merge matches across both strands
        if (!strand_failed && !matches.empty()) {
          matches_by_strand.insert(matches_by_strand.end(), 
            matches.begin(), matches.end());
        }

        // '+' strand: continue to try the other strand
        if (strand == '+') continue;

        // '-' strand: finished both, return combined matches
        return {true, matches_by_strand};
      }
    }

    // If this was an unstranded read and both strands fail, 
    // give up
    return config.default_result;
  }

  ReadEvaluationResult 
  ShortReadEvaluator::evaluate(BundleData *bundle, read_id_t id, g2tTree *g2t) {

    uint32_t boundary_tolerance = 5; 
    uint32_t max_clip = 25;
    uint32_t max_ins = 0;
    uint32_t max_gap = 0;
    double similarity_threshold = 0.90;

    ReadEvaluationConfig config = {
      boundary_tolerance,             // boundary tolerance
      max_clip,                       // max clip size
      max_ins,                        // max insertion to intervals
        // for when there exist no guides to explain a portion of an exon
      max_gap,                        // max gap in reference to intervals
      true,                           // filter similarity?
      similarity_threshold,           // similarity threshold
      false,                          // ignore small exons?
      0,                              // small exon size
      {false, {}}                     // default result
    }; 

    ReadEvaluationResult res = evaluate_exon_chains(bundle, id, g2t, config);
    return res;
  }

  ReadEvaluationResult 
  LongReadEvaluator::evaluate(BundleData *bundle, read_id_t id, g2tTree *g2t) {

    uint32_t boundary_tolerance = 20;
    uint32_t max_clip = 20;
    uint32_t max_ins = 20;
    uint32_t max_gap = 20;
    double similarity_threshold = 0.99;
    ReadEvaluationConfig config = {
      boundary_tolerance,             // boundary tolerance
      max_clip,                       // max clip size
      max_ins,                        // max insertion to intervals
        // for when there exist no guides to explain a portion of an exon
      max_gap,                        // max gap in reference to intervals
      true,                           // filter similarity?
      similarity_threshold,           // similarity threshold
      true,                           // ignore small exons?
      SMALL_EXON,                     // small exon size
      {false, {}}                     // default result
    };   

    ReadEvaluationResult res = evaluate_exon_chains(bundle, id, g2t, config); 
    return res;
  }

} // namespace bramble