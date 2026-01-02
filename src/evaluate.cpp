
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
  std::unordered_map<tid_t, ExonChainMatch> 
  ReadEvaluator::evaluate(CReadAln * read, 
                          read_id_t id, g2tTree *g2t) {
    return std::unordered_map<tid_t, ExonChainMatch>{};
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

  // std::tuple<uint32_t, uint32_t> 
  // ReadEvaluator::get_exon_coordinates(GSeg exon, char strand, 
  //                                     uint32_t ref_len) {
  //   uint32_t exon_start, exon_end;
  //   //if (strand == '+') {
  //     exon_start = exon.start;
  //     exon_end = exon.end + 1;
  //   // } else {
  //   //   exon_start = ref_len - (exon.end + 1);
  //   //   exon_end = ref_len - exon.start;
  //   // }
  //   return std::make_tuple(exon_start, exon_end);
  // }
    
  // /**
  // * Get the match position for this transcript
  // *
  // * @param interval first interval the read matched to
  // * @param tid transcript id
  // * @param read_strand strand that the read aligned to
  // * @param g2t g2t tree for bundle
  // * @param exon_start start of first read exon
  // * (USE_FASTA mode)
  // */
  // pos_t ReadEvaluator::get_match_pos(std::shared_ptr<IntervalNode> interval, char strand,
  //                                   tid_t tid, uint32_t exon_start, uint32_t exon_end) {
    
  //   pos_t cum_len = interval->tid_cum_len[tid];
  //   uint32_t match_start;
  //   if (strand == '+') {
  //     match_start = (exon_start > interval->start) ? 
  //       (exon_start - interval->start) : 0;
  //   } else {
  //     match_start = (exon_end < interval->end) ? 
  //       (interval->end - exon_end) : 0;
  //   }
  //   return (pos_t)(match_start + cum_len);
  // }

  std::vector<char> 
  ReadEvaluator::get_strands_to_check(CReadAln * read) {
    if (read->strand == '+') return {'+'}; 
    if (read->strand == '-') return {'-'}; 
    return {'+', '-'};  // unstranded: try forward first, then reverse
  }

  // Add operation to every CIGAR in matches
  // Used for hide_small_exon
  void ReadEvaluator::add_operation_for_all(std::unordered_map<tid_t, ExonChainMatch> &matches,
                                            uint32_t length, uint8_t op, g2tTree *g2t, 
                                            char strand) {
    auto it = matches.begin();
    while (it != matches.end()) {
      std::shared_ptr<Cigar> match_cigar = it->second.align.cigar;
      match_cigar->add_operation(length, op);
      if (op != BAM_CSOFT_CLIP && op != BAM_CINS) {
        it->second.ref_consumed += length;
        it->second.total_coverage += length;
      }
      ++it;
    }
  };

  void ReadEvaluator::hide_small_exon(std::unordered_map<tid_t, ExonChainMatch> &matches,
                                      uint32_t exon_start, uint32_t exon_end,
                                      ExonStatus status, g2tTree *g2t, char strand) {

    uint32_t match_size = exon_end - exon_start;

    // First exon: soft clip whole thing
    if ((status == FIRST_EXON || status == ONLY_EXON)) {
      add_operation_for_all(matches, match_size, BAM_CSOFT_CLIP, g2t, strand);
    }
    
    // Middle exon: entire thing is an insertion
    if ((status == MIDDLE_EXON)) {
      add_operation_for_all(matches, match_size, BAM_CINS, g2t, strand);
    }
    
    // Last exon: soft clip whole thing
    if ((status == LAST_EXON || status == ONLY_EXON)) {
      add_operation_for_all(matches, match_size, BAM_CSOFT_CLIP, g2t, strand);
    }
  }
  
  // goal: remove tids with intermediate exon unaccounted for
  void 
  ReadEvaluator::ensure_continuity(std::vector<std::shared_ptr<IntervalNode>> &intervals,
                                  std::unordered_map<tid_t, uint8_t> &exon_id_map,
                                  ExonStatus status, char strand) {
    
    // Keep only TIDs with correct continuity
    if (status == MIDDLE_EXON || status == LAST_EXON) {
      auto it_i = intervals.begin();
      while (it_i != intervals.end()) {
        auto interval = *it_i;
        auto tid = interval->tid;
        
        // Check continuity: successive query exons should match to 
        // successive guide exons, with no guide exons in between

        if (!exon_id_map.count(tid)) {
          it_i = intervals.erase(it_i);
          continue;
        }

        uint8_t prev_exon_id = exon_id_map[tid];
        uint8_t curr_exon_id = interval->exon_id;

        if (strand == '+' && curr_exon_id != (prev_exon_id + 1)) {
          it_i = intervals.erase(it_i);
        } else if (strand == '-' && curr_exon_id != (prev_exon_id - 1)) {
          it_i = intervals.erase(it_i);
        } else {
          ++it_i;
        }
      }
    }

    exon_id_map.clear();
    for (auto &interval : intervals) {
      auto tid = interval->tid;
      exon_id_map[tid] = interval->exon_id;
    }
  }

  void ReadEvaluator::build_exon_cigar(std::unordered_map<tid_t, ExonChainMatch> &matches,
                                      ExonStatus status, uint32_t exon_start, uint32_t exon_end,
                                      std::vector<std::shared_ptr<IntervalNode>> &intervals,
                                      std::unordered_map<tid_t, uint32_t> &gaps,
                                      g2tTree* g2t, char strand) {
    // Handle gaps between previous and current exon
    if (!gaps.empty()) {
      for (auto &pair : matches) {
        tid_t tid = pair.first;
        auto match = pair.second;
        std::shared_ptr<Cigar> cigar = match.align.cigar;
        if (gaps.count(tid) && gaps[tid] > 0) cigar->add_operation(gaps[tid], BAM_CDEL);
        match.ref_consumed += gaps[tid];
        match.total_coverage += gaps[tid];
      }
    }

    // Loop through intervals
    // Add S, M, I to match CIGAR accordingly
    auto it_i = intervals.begin();
    while (it_i != intervals.end()) {
      auto interval = *it_i;
      auto tid = interval->tid;
      std::string tid_string = g2t->getTidName(tid);
      // if (tid_string == "CHS.55407.18" && exon_start == 35682013) {
      //   printf("interval->start = %d, interval->end = %d\n", interval->start, interval->end);
      //   printf("exon_start = %d, exon_end = %d\n", exon_start, exon_end);
      // }
      if (!matches.count(tid)) {
        it_i = intervals.erase(it_i);
        continue;
      }
      auto& match = matches[tid];
      std::shared_ptr<Cigar> cigar = match.align.cigar;

      if (strand == '+') {
        // Handle start boundary
        if (exon_start < interval->start) {
          uint32_t overhang = interval->start - exon_start;
          // if (tid_string == "CHS.55407.18" && exon_start == 35682013) {
          //   printf("left clip = %d\n", overhang);
          // }
          if (status == FIRST_EXON || status == ONLY_EXON) 
            cigar->add_operation(overhang, BAM_CSOFT_CLIP);
            // if (tid_string == "CHS.55407.18" && exon_start == 35682013) {
            //   printf("added clip\n");
            // }
          else if (status == MIDDLE_EXON || status == LAST_EXON)
            cigar->add_operation(overhang, BAM_CINS);
        }
      } else {
        // Handle end boundary
        if (exon_end > interval->end) {
          uint32_t overhang = exon_end - interval->end;
          if (status == LAST_EXON || status == ONLY_EXON) 
            cigar->add_operation(overhang, BAM_CSOFT_CLIP);
          else if (status == FIRST_EXON || status == MIDDLE_EXON)
            cigar->add_operation(overhang, BAM_CINS);
        }
      }
 
      // Add match for overlapping region
      uint32_t overlap_start = std::max(exon_start, interval->start);
      uint32_t overlap_end = std::min(exon_end, interval->end);  
      if (overlap_end >= overlap_start) {
        uint32_t match_length = overlap_end - overlap_start;
        // if (tid_string == "CHS.55407.18" && exon_start == 35682013) {
        //   printf("match length = %d\n", match_length);
        // }
        cigar->add_operation(match_length, BAM_CMATCH);
        match.ref_consumed += match_length;
        match.total_coverage += match_length;
      }

      if (strand == '+') {
        // Handle end boundary
        if (exon_end > interval->end) {
          uint32_t overhang = exon_end - interval->end;
          if (status == LAST_EXON || status == ONLY_EXON) 
            cigar->add_operation(overhang, BAM_CSOFT_CLIP);
          else if (status == FIRST_EXON || status == MIDDLE_EXON)
            cigar->add_operation(overhang, BAM_CINS);
        }
      } else {
        // Handle start boundary
        if (exon_start < interval->start) {
          uint32_t overhang = interval->start - exon_start;
          // if (tid_string == "CHS.55407.18" && exon_start == 35682013) {
          //   printf("left clip = %d\n", overhang);
          // }
          if (status == FIRST_EXON || status == ONLY_EXON) 
            cigar->add_operation(overhang, BAM_CSOFT_CLIP);
            // if (tid_string == "CHS.55407.18" && exon_start == 35682013) {
            //   printf("added clip\n");
            // }
          else if (status == MIDDLE_EXON || status == LAST_EXON)
            cigar->add_operation(overhang, BAM_CINS);
        }
      }

      ++it_i;
    }
  }

  void ReadEvaluator::compile_matches(std::unordered_map<tid_t, ExonChainMatch> &matches,
                                      std::vector<std::shared_ptr<IntervalNode>> &intervals,
                                      char strand, uint32_t exon_start, uint32_t exon_end,
                                      uint32_t cumulative_length, uint32_t coverage_augment) {
    // Initialize matches with TIDs and positions
    for (const auto &interval : intervals) {
      ExonChainMatch match;
      match.align.pos = interval->pos;
      match.align.strand = strand;
      match.align.cigar = std::make_shared<Cigar>();
      match.ref_consumed = 0;
      match.total_coverage = 0 + coverage_augment; // same as ref consumed
      match.total_read_length = cumulative_length + (exon_end - exon_start);
      match.align.similarity_score = 0;
      matches[interval->tid] = match;

      if (cumulative_length > 0) {
        match.align.cigar->add_operation(cumulative_length, BAM_CSOFT_CLIP);
      }
    }
  }

  bool ReadEvaluator::update_matches(std::unordered_map<tid_t, ExonChainMatch> &matches,
                                    std::vector<std::shared_ptr<IntervalNode>> &intervals, 
                                    char strand, uint32_t exon_start, uint32_t exon_end, 
                                    ExonStatus status) {
    
    std::unordered_map<tid_t, pos_t> candidate_tids;
    for (auto &interval : intervals) {
      candidate_tids[interval->tid] = interval->pos;
    }

    auto it_m = matches.begin();
    while (it_m != matches.end()) {
      auto match = &(it_m->second);
      match->total_read_length += (exon_end - exon_start);

      if (status == MIDDLE_EXON || status == LAST_EXON) {
        tid_t tid = (*it_m).first;

        // Only keep tids supported by all exons
        if (candidate_tids.count(tid)) {

          // For reads on negative strand, need to update
          // match position at each new exon
          if (strand == '-') {
            pos_t pos = candidate_tids[tid];
            match->align.pos = pos;
          }
          ++it_m;

        } else {
          it_m = matches.erase(it_m);
        }
 
      } else {
        ++it_m;
      }
      
    }
    return !matches.empty();
  }

  int32_t get_rand(uint32_t x) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, x - 1);
    return dis(gen);
  }

  void 
  ReadEvaluator::filter_by_similarity(std::unordered_map<tid_t, ExonChainMatch> &matches, 
                                      double similarity_threshold) {

    for (auto it = matches.begin(); it != matches.end(); ) {
      auto& match = it->second;

      double similarity = (match.total_read_length > 0)
        ? (static_cast<double>(match.total_coverage) / match.total_read_length)
        : 0.0;

      if (similarity > similarity_threshold) {
        match.align.similarity_score = similarity;
        ++it;
      } else {
        it = matches.erase(it);
      }
    }

    if (matches.empty()) return;

    auto best_it = matches.end();
    double best_score = -1.0;
    int32_t hit_index = 1;

    for (auto it = matches.begin(); it != matches.end(); ++it) {
      auto& match = it->second;
      match.align.hit_index = hit_index++;

      if (match.align.similarity_score > best_score) {
        best_score = match.align.similarity_score;
        best_it = it;
      }
    }

    if (best_it != matches.end()) {
      best_it->second.align.primary_alignment = true;
    } else {
      int32_t rand_idx = get_rand(matches.size());
      auto it = matches.begin();
      std::advance(it, rand_idx);
      it->second.align.primary_alignment = true;
    }
  }

  std::unordered_map<tid_t, ExonChainMatch> 
  ReadEvaluator::evaluate_exon_chains(CReadAln *read, 
                                      read_id_t id, g2tTree *g2t, 
                                      ReadEvaluationConfig config) {
    GVec<GSeg> read_exons = read->segs;
    uint32_t exon_count = read_exons.Count();
    refid_t refid = read->refid;
    std::string read_name = read->brec->name();

    std::unordered_map<tid_t, ExonChainMatch> matches;
    std::unordered_map<tid_t, ExonChainMatch> matches_by_strand;
    //std::vector<std::shared_ptr<IntervalNode>> intervals;
    std::vector<uint32_t> exon_sizes;
    std::unordered_map<tid_t, uint8_t> exon_id_map;

    // printf("read_name = %s\n", read_name.c_str());

    auto strands_to_check = get_strands_to_check(read);
    for (char strand : strands_to_check) {
      matches.clear();
      //intervals.clear();
      exon_sizes.clear();
      exon_id_map.clear();
      bool strand_failed = false;
      uint32_t total_exon_matches = exon_count;
      uint32_t cumulative_length = 0;
      uint32_t coverage_augment = 0;

      for (uint j = 0; j < exon_count; j++) {
        auto exon = read_exons[j];
        exon.end++; // make exclusive

        // if (read_name == "read85009325/CHS.44955.2;mate1:613-713;mate2:613-713") {
        //   printf("exon.start = %d, exon.end = %d\n", exon.start, exon.end);
        //   config.print = true;
        // } else {
        //   config.print = false;
        // }

        ExonStatus status;
        if (exon_count == 1) {
          status = ONLY_EXON;
        } else if (j == 0) {
          status = FIRST_EXON;
        } else if (j < (exon_count - 1)) {
          status = MIDDLE_EXON;
        } else {
          status = LAST_EXON;  
        }

        bool is_small_exon = (exon.end - exon.start <= config.small_exon_size);

        // Find all guide intervals that contain the read exon
        auto intervals = g2t->getIntervals(refid, strand, exon,
          config, status);
        if (intervals.empty()) {
          // printf("INTERVALS EMPTY\n");
          // printf("read_name == %s\n", read_name.c_str());
          // printf("exon_i = %d\n", j);
          // printf("refid = %d\n", refid);
          // printf("exon.start = %d, exon.end = %d\n", exon.start, exon.end);
          // printf("\n");
          if (config.ignore_small_exons & is_small_exon) {
            hide_small_exon(matches, exon.start, exon.end, 
              status, g2t, strand); 
            total_exon_matches--;
            cumulative_length += (exon.end - exon.start);
            //coverage_augment += (exon.end - exon.start);
            continue;
          }
          dropped_reads++;
          strand_failed = true;
          matches.clear();
          break;
        }

        // Filter tids
        ensure_continuity(intervals, exon_id_map, status, strand);
        if (intervals.empty()) {
          strand_failed = true;
          matches.clear();
          if (config.print) {
            printf("Strand failed at filter tids\n");
          }
          
          unresolved_reads++;
          break;
        }

        if (config.print) {
          for (auto &pair : exon_id_map) {
            tid_t tid = pair.first;
            uint8_t ei = pair.second;
            printf("tid = %d, ei = %d\n", tid, ei);
          }
        }

        // Create exon chain matches 
        if (status == FIRST_EXON || status == ONLY_EXON) {
          compile_matches(matches, intervals, 
            strand, exon.start, exon.end, cumulative_length, 
            coverage_augment);
        }

        // Build CIGAR as we process each exon
        std::unordered_map<tid_t, uint32_t> gaps;
        build_exon_cigar(matches, status, exon.start, exon.end, 
          intervals, gaps, g2t, strand);

        if (status == MIDDLE_EXON || status == LAST_EXON) {
          strand_failed = !update_matches(matches, intervals, strand, 
          exon.start, exon.end, status);
          if (strand_failed) {
            matches.clear();
            if (config.print) { 
              printf("Strand failed at update matches\n");
            }
            
            unresolved_reads++;
            break;
          }
        }

        // Fail if we run out of candidate TIDs
        if (intervals.empty()) {
          unresolved_reads++;
          if (config.print) {
            printf("Strand failed at candidate tids empty\n");
          }
          
          strand_failed = true;
          matches.clear();
          break;
        }

        // After the last exon, filter matches by similarity
        if (config.filter_similarity && 
          (status == ONLY_EXON || status == LAST_EXON)) 
          filter_by_similarity(matches, config.similarity_threshold);

        // if (read_name == "read3599125/CHS.2231.5;mate1:288-388;mate2:288-388") {
          //printf("tids after filter similarity \n");
        // for (auto &pair : matches) {
        //   std::string tid_string = g2t->getTidName(pair.first);
        //   if (tid_string == "CHS.55407.18" && exon.start == 35682013) {
        //     printf("read name = %s\n", read_name.c_str());
        //     printf("tid = %d, %s\n", pair.first, tid_string.c_str());
        //     printf("pos = %d\n", pair.second.align.pos);
        //     printf("IDEAL CIGAR: ");
        //     for (const auto& pair : pair.second.align.cigar->cigar) {
        //       printf("%u%c ", pair.first, "MIDNSHP=XB"[pair.second]);
        //     }
        //     printf("\n\n");
        //   }
            
        //   // }
        // }
        
        if ((status == ONLY_EXON || status == LAST_EXON) 
          && matches.empty()) {
          unresolved_reads++;
          if (config.print) {
            printf("Strand failed at filter similarity\n");
          }
          
          strand_failed = true;
          matches.clear();
          break;
        }

        cumulative_length += (exon.end - exon.start);
        if (total_exon_matches == 0) {
          unresolved_reads++;
          strand_failed = true;
          matches.clear();
          break;
        }

      } // end of exon loop

      // Handle short reads
      if (!LONG_READS) {
        // Worked on this strand, so we are done
        if (!strand_failed && !matches.empty()) {
          return matches;
        }
        // Stranded read failed -> fail
        if (read->strand != '.') {
          return matches;
        }
        // Unstranded: try the opposite strand next
        continue;
      }

      // Handle long reads
      if (LONG_READS) {
        // Stranded long reads: either success or failure
        if (read->strand != '.') {
          if (!strand_failed && !matches.empty()) {
            return matches;
          } else {
            return matches;
          }
        }

        // Unstranded long reads: merge matches across both strands
        if (!strand_failed && !matches.empty()) {
          for (auto &pair : matches) {
            tid_t tid = pair.first;
            auto match = pair.second;
            matches_by_strand[tid] = match;
          }
        }

        // '+' strand: continue to try the other strand
        if (strand == '+') continue;

        // '-' strand: finished both, return combined matches
        return matches_by_strand;
      }
    }

    // If this was an unstranded read and both strands fail, 
    // give up
    matches.clear();
    return matches;
  }

  std::unordered_map<tid_t, ExonChainMatch> 
  ShortReadEvaluator::evaluate(CReadAln * read, 
                              read_id_t id, g2tTree *g2t) {

    uint32_t boundary_tolerance = 5; //5
    uint32_t max_clip = 5; //25
    uint32_t max_ins = 5;
    uint32_t max_gap = 5;
    uint32_t max_junc_gap = 0;
    double similarity_threshold = 0.90; // 0.90

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
      max_junc_gap,                   // max junction gap
      false                           // print debug statements
    }; 

    auto matches = evaluate_exon_chains(read, id, g2t, config);
    // for (auto &pair : matches) {
    //   tid_t tid = pair.first;
    //   auto match = pair.second;
    //   std::string tid_string = g2t->getTidName(tid);
    //   printf("match tid = %s\n", tid_string.c_str());
    // }
    return matches;
  }

  std::unordered_map<tid_t, ExonChainMatch> 
  LongReadEvaluator::evaluate(CReadAln * read, 
                              read_id_t id, g2tTree *g2t) {

    uint32_t boundary_tolerance = 20;
    uint32_t max_clip = 20;
    uint32_t max_ins = 20;
    uint32_t max_gap = 20;
    uint32_t max_junc_gap = 0;
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
      max_junc_gap,                   // max junction gap
      false                           // print debug statements
    };   

    auto matches = evaluate_exon_chains(read, id, g2t, config); 
    return matches;
  }

} // namespace bramble