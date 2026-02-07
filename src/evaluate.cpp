
#include <algorithm>
#include <memory>
#include <set>
#include <string>
#include <vector>
#include <unordered_set>
#include <cmath>
#include <random>
#include <list>

#include "types.h"
#include "bramble.h"
#include "g2t.h"
#include "evaluate.h"
#include "bam.h"
#include "sw.h"

#ifndef NOTHREADS
#include "GThreads.h"
#endif

extern GFastMutex bam_io_mutex;  // protects BAM io

extern bool LONG_READS;
extern bool USE_FASTA;
extern bool SOFT_CLIPS;
extern bool STRICT;

extern uint32_t dropped_reads;
extern uint32_t unresolved_reads;

namespace bramble {

  ReadEvaluator::~ReadEvaluator() {}

  // Evaluate read group and return matches
  std::unordered_map<tid_t, ExonChainMatch> 
  ReadEvaluator::evaluate(CReadAln * read, 
                          read_id_t id, 
                          std::shared_ptr<g2tTree> g2t,
                          uint8_t *seq, int seq_len) {
    return std::unordered_map<tid_t, ExonChainMatch>{};
  }

  ExonStatus get_exon_status(uint32_t exon_count, uint32_t j) {
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
    return status;
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

  std::vector<char> 
  ReadEvaluator::get_strands_to_check(CReadAln *read) {
    if (LONG_READS && read->segs.Count() == 1) {
      return {'+', '-'};  // correct for bad strand guess
    }
    if (read->strand == '+') return {'+'}; 
    if (read->strand == '-') return {'-'}; 
    return {'+', '-'};  // unstranded: try forward first, then reverse
  }

  void ReadEvaluator::get_clips(CReadAln *read, ReadEvaluationConfig &config,
                                bool &failure, bool &has_left_clip, bool &has_right_clip,
                                uint32_t &n_left_clip, uint32_t &n_right_clip) {
    auto brec = read->brec;
    bam1_t* b = read->brec->get_b();
    uint32_t* cigar = bam_get_cigar(b);
    uint32_t n_cigar = b->core.n_cigar;

    if (!cigar) {
      printf("Error: alignment is missing CIGAR string\n");
      printf("Skipping this alignment\n");
      failure = true;
      return;
    }
    
    if ((cigar[0] & BAM_CIGAR_MASK) == BAM_CHARD_CLIP) {
      if (1 < n_cigar && cigar[1] & BAM_CIGAR_MASK == BAM_CSOFT_CLIP) {
        has_left_clip = USE_FASTA;
        n_left_clip = cigar[1] >> BAM_CIGAR_SHIFT;
      }
       
    } else if ((cigar[0] & BAM_CIGAR_MASK) == BAM_CSOFT_CLIP) {
      has_left_clip = USE_FASTA;
      n_left_clip = cigar[0] >> BAM_CIGAR_SHIFT;
    }
      
    if ((cigar[n_cigar - 1] & BAM_CIGAR_MASK) == BAM_CHARD_CLIP) {
      if (n_cigar - 2 >= 0 && cigar[n_cigar - 2] & BAM_CIGAR_MASK == BAM_CSOFT_CLIP) {
        has_right_clip = USE_FASTA;
        n_right_clip = cigar[n_cigar - 2] >> BAM_CIGAR_SHIFT;
      }
        
    } else if ((cigar[n_cigar - 1] & BAM_CIGAR_MASK) == BAM_CSOFT_CLIP) {
      has_right_clip = USE_FASTA;
      n_right_clip = (cigar[n_cigar - 1] >> BAM_CIGAR_SHIFT);
    }

    if (config.print) {
      printf("Has left clip? %d # left clip = %d\n", has_left_clip, n_left_clip);
      printf("Has right clip? %d # right clip = %d\n", has_right_clip, n_right_clip);
    }
  }



  void ReadEvaluator::get_intervals(std::unordered_map<tid_t, TidData> &data,
                                    CReadAln *read, uint32_t j, uint32_t exon_count,
                                    ReadEvaluationConfig &config, std::shared_ptr<g2tTree> g2t,
                                    refid_t refid, char strand, bool has_left_clip,
                                    bool has_right_clip, Segment &start_ignore,
                                    Segment &end_ignore, bool &has_gap, bool &failure) {
    auto qexon = read->segs[j];

    auto status = get_exon_status(exon_count, j); 
      // first, middle, last, or only exon?

    bool is_small_exon = (qexon.end - qexon.start <= config.small_exon_size);
    bool data_empty = data.empty();

    // Find all guide exons that contain the query exon
    auto guide_exons = g2t->getGuideExons(refid, strand, 
      qexon, config, status);

    std::set<tid_t> candidate_tids;
    if (!guide_exons.empty()) {

      int prev_start = -1;
      int prev_end = -1;

      for (auto &pair : guide_exons) {
        tid_t tid = pair.first;
        auto &gexon = pair.second;

        if (prev_start == -1) {
          prev_start = gexon->start;
          prev_end = gexon->end;
        } else if (gexon->start != prev_start
          || gexon->end != prev_end) {
          //has_variants = true;
        }

        candidate_tids.emplace(tid);

        if (data_empty && status != INS_EXON) {
          TidData tid_data;
          tid_data.has_left_clip = has_left_clip;
          tid_data.has_right_clip = has_right_clip;
          tid_data.has_left_ins = (gexon->left_ins > 0);
          tid_data.n_left_ins = gexon->left_ins;
          if (config.print) printf("gexon->left_ins = %d\n", gexon->left_ins);
          data[tid] = tid_data;

          if (gexon->left_gap > 0) has_gap = true;
        } else {
          if (!data.count(tid) || data[tid].elim) continue;
        }

        if (status == LAST_EXON || status == ONLY_EXON) {
          auto &tid_data = data[tid];
          tid_data.has_right_ins = (gexon->right_ins > 0);
          tid_data.n_right_ins = gexon->right_ins;
          if (config.print) printf("gexon->right_ins = %d\n", gexon->right_ins);

          if (gexon->right_gap > 0) has_gap = true;
        }

        Segment segment;
        segment.is_valid = true;
        segment.qexon = qexon;
        segment.gexon = gexon;
        segment.has_qexon = true;
        segment.has_gexon = true;
        segment.status = status;
        segment.is_small_exon = is_small_exon;
        data[tid].segments.emplace_back(segment);
      }

      for (auto &pair : data) {
        tid_t tid = pair.first;
        if (!candidate_tids.count(tid)) {
          pair.second.elim = true;
        }
      }

      return;
    }

    if (config.print) printf("\nGuide exons empty\n");
    if (status != ONLY_EXON && config.ignore_small_exons 
      && is_small_exon) {

      // this whole block technically allows any instance of 'guide exons empty'
      // even if there were guide exons and they were just bad
    
      if (status == MIDDLE_EXON) {
        if (data.empty()) {
          dropped_reads++;
          failure = true;
          return;
        }
        for (auto &pair: data) {
          tid_t tid = pair.first;
          auto &tid_data = pair.second;
          
          Segment ignore;
          ignore.is_valid = true;
          ignore.qexon = qexon;
          ignore.has_qexon = true;
          ignore.has_gexon = false;
          ignore.status = INS_EXON;
          ignore.is_small_exon = true;
          tid_data.segments.emplace_back(ignore);
        }

        // these didn't help
      // } else if (status == FIRST_EXON) {
      //   start_ignore.is_valid = true;
      //   start_ignore.qexon = qexon;
      //   start_ignore.has_qexon = true;
      //   start_ignore.has_gexon = false;
      //   start_ignore.status = INS_EXON;
      //   start_ignore.is_small_exon = true;
      // } else if (status == LAST_EXON) {
      //   end_ignore.is_valid = true;
      //   end_ignore.qexon = qexon;
      //   end_ignore.has_qexon = true;
      //   end_ignore.has_gexon = false;
      //   end_ignore.status = INS_EXON;
      //   end_ignore.is_small_exon = true;
      // }
      } else {
        failure = true;
        return;
      }

      return; // move on to next query exon
    }

    // no resolution; give up on read
    failure = true;
    return;
  }

  void ReadEvaluator::correct_for_gaps(TidData &tid_data, tid_t tid,
                                      ReadEvaluationConfig &config, 
                                      std::shared_ptr<g2tTree> g2t,
                                      char strand, refid_t refid) {
    auto it_s = tid_data.segments.begin();
    std::shared_ptr<GuideExon> prev_gexon = nullptr;
    
    bool first = true;
    while (it_s != tid_data.segments.end()) {
      auto &segment = *it_s;

      if (!segment.has_gexon) {
        ++it_s;
        continue;
      }

      auto gexon = segment.gexon;
      
      // Skip continuity check for first exon
      if (first) {
        first = false;
        prev_gexon = gexon;
        ++it_s;
        continue;
      }

      uint8_t gap = gexon->exon_id - prev_gexon->exon_id;

      // Short reads: require strict continuity
      if (!LONG_READS) {
        if (gap != 1) {
          tid_data.elim = true;
          return;
        }
        prev_gexon = gexon;
        ++it_s;
        continue;
      }

      // Long reads: handle gaps
      if (gap > 2) {
        tid_data.elim = true;
        return;
      }
      
      if (gap == 2) {
        bool is_plus = (strand == '+');
        uint32_t gap_start = is_plus ? gexon->prev_start : gexon->next_start;
        uint32_t gap_end = is_plus ? gexon->prev_end : gexon->next_end;

        if ((gap_start == 0 && gap_end == 0) || 
          (gap_end - gap_start > config.max_gap)) {
          tid_data.elim = true;
          return;
        }

        auto prev_exon = g2t->getGuideExonForTid(refid, strand, tid, 
          gap_start, gap_end);
        if (!prev_exon) {
          tid_data.elim = true;
          return;
        }

        Segment gap_seg;
        gap_seg.is_valid = true;
        gap_seg.gexon = prev_exon;
        gap_seg.has_gexon = true;
        gap_seg.status = GAP_EXON;
        gap_seg.is_small_exon = (prev_exon->end - prev_exon->start 
          <= config.small_exon_size);
        
        it_s = tid_data.segments.emplace(it_s, gap_seg);
        ++it_s;

        if (config.print) {
            printf("Inserted GAP: tid=%d exon_id=%d [%d-%d]\n",
                  tid, prev_exon->exon_id, prev_exon->start, prev_exon->end);
        }
      }
      
      prev_gexon = gexon;
      ++it_s;
    }
  }

  void ReadEvaluator::left_clip_rescue(TidData &tid_data, Segment &start_ignore,
                                        char strand, std::shared_ptr<g2tTree> g2t,
                                        refid_t refid, tid_t tid, uint32_t n_left_clip,
                                        uint32_t n_left_ins, ReadEvaluationConfig &config, 
                                        CReadAln *read, uint8_t *seq, int seq_len) {
    auto it = tid_data.segments.begin();
    if (start_ignore.is_valid) {
      tid_data.has_left_clip = false;
      return;
    }
    Segment &segment = *it;

    auto gexon = segment.gexon;
    if (!segment.has_gexon) {
      tid_data.has_left_clip = false;
      return;
    }

    if (config.print) printf("gexon->left_gap = %d\n", gexon->left_gap);
    if (gexon->left_gap > 0) {
      tid_data.has_left_clip = false;
      return;
    }

    bool use_same_exon = false;
    // if (gexon->left_gap >= 0) {
    //   use_same_exon = true;
    // }

    int start_pos = 0;  // Starting position of soft clip
    uint32_t total_bases = n_left_clip + n_left_ins;

    // Extract the soft-clipped bases
    char rescue_seq[total_bases + 1];
    for (int i = 0; i < total_bases; i++) {
      rescue_seq[i] = seq_nt16_str[bam_seqi(seq, start_pos + i)];
    }
    rescue_seq[total_bases] = '\0';  // Null terminate

    std::string remaining_qseq = rescue_seq;
    std::shared_ptr<GuideExon> current_gexon = gexon;
    
    // Safety limit to prevent infinite loops
    const int MAX_EXONS = 20;
    int exons_processed = 0;

    while (!remaining_qseq.empty() && exons_processed < MAX_EXONS) {
      exons_processed++;
      
      // Get previous guide exon (moving leftward)
      std::shared_ptr<GuideExon> left_gexon;
      std::string next_gseq;
      if (use_same_exon && exons_processed == 1) {
        left_gexon = current_gexon;
        int size = (segment.qexon.start - left_gexon->start);
        next_gseq = (left_gexon->seq).substr(0, size);
      }
      else if (strand == '+') {
        if (!current_gexon->has_prev) {
          if (exons_processed == 1) {
            tid_data.has_left_clip = false;
            return;
          } else {
            break;
          }
        }
        left_gexon = g2t->getGuideExonForTid(refid, strand, tid, 
          current_gexon->prev_start, current_gexon->prev_end);
        next_gseq = left_gexon->seq;
      } else {
        if (!current_gexon->has_next) {
          if (exons_processed == 1) {
            tid_data.has_left_clip = false;
            return;
          } else {
            break;
          }
        }
        left_gexon = g2t->getGuideExonForTid(refid, strand, tid, 
          current_gexon->next_start, current_gexon->next_end);
        next_gseq = left_gexon->seq;
      }

      if (config.print) {
        auto tid_string = g2t->getTidName(tid);
        printf("\nGetting left alignment for tid = %s (exon %d)\n", 
               tid_string.c_str(), exons_processed);
        printf("remaining_qseq = %s (len=%zu)\n", remaining_qseq.c_str(), remaining_qseq.length());
        printf("next_gseq = %s (len=%zu)\n", next_gseq.c_str(), next_gseq.length());
        printf("curr_gexon_start = %d, curr_gexon_end = %d\n", current_gexon->start, current_gexon->end);
        printf("next_gexon_start = %d, next_gexon_end = %d\n", left_gexon->start, left_gexon->end);
      }

      AlignmentResult result = smith_waterman(remaining_qseq, next_gseq, END);
        // returns INCLUSIVE end values

      if (config.print) {
        printf("\nLEFT CLIP ALIGNMENT RESULT:\n");
        printf("score = %d\n", result.score);
        printf("result.accuracy = %f\n", result.accuracy);
        printf("query_start = %d, query_end = %d\n", result.start_i, result.end_i);
        printf("ref_start = %d, ref_end = %d\n", result.start_j, result.end_j);
      }

      if (exons_processed == 1 && (result.accuracy <= 0.70 
        || result.score < 10)) {
        tid_data.has_left_clip = false;
        return;
      }

      // Query exon coordinates based on where alignment occurs in guide
      uint32_t qstart = left_gexon->start + result.start_j;
      uint32_t qend = left_gexon->start + result.end_j;
      
      left_gexon->pos = result.pos + left_gexon->pos_start;

      if (exons_processed == 1) {
        // If we have an insertion, correct the current segment
        if (n_left_ins > 0) {
          segment.qexon.start = gexon->start; // remove insertion from this exon
        }
        if (use_same_exon) { 
          segment.qexon.start = gexon->start; // remove deletion from this exon
        }
      }

      Segment left_clip;
      left_clip.is_valid = true;
      left_clip.qexon = GSeg(qstart, qend);
      left_clip.gexon = left_gexon;
      left_clip.has_qexon = true;
      left_clip.has_gexon = true;
      left_clip.status = LEFTC_EXON;
      left_clip.is_small_exon = (remaining_qseq.length() <= config.small_exon_size);
      left_clip.cigar = result.cigar;
      left_clip.score = result.score;
      tid_data.segments.emplace(tid_data.segments.begin(), left_clip);
      tid_data.has_left_clip = true;

      // Calculate remaining query sequence - moving leftward
      // For left rescue: alignment consumed from start_i to end_i (rightmost part)
      // We keep everything before start_i (leftmost part) for next iteration
      
      // Always update remaining_qseq first
      if (result.start_i > 1) {
        // should add deletion here

        remaining_qseq = remaining_qseq.substr(0, result.start_i);
      } else {
        remaining_qseq.clear();  // Everything consumed
      }

      // Then check if we should continue
      if (remaining_qseq.empty()) {
        break;
      }

      // Update for next iteration
      current_gexon = left_gexon;
      break;
    }

    // After while loop ends
    if (!remaining_qseq.empty()) {
      // Add final soft clip for any leftover bases
      tid_data.segments.front().cigar.prepend_operation(remaining_qseq.length(), BAM_CLIP_OVERRIDE);
    }
  }

  void ReadEvaluator::right_clip_rescue(TidData &tid_data, Segment &end_ignore,
                                        char strand, std::shared_ptr<g2tTree> g2t,
                                        refid_t refid, tid_t tid, uint32_t n_right_clip,
                                        uint32_t n_right_ins, ReadEvaluationConfig &config, 
                                        CReadAln *read, uint8_t *seq, int seq_len) {
    auto rit = tid_data.segments.rbegin(); // last element
    if (end_ignore.is_valid) {
      tid_data.has_right_clip = false;
      return;
    }
    Segment &segment = *rit;

    auto gexon = segment.gexon;
    if (!segment.has_gexon) {
      if (config.print) printf("No exon in segment\n");
      tid_data.has_right_clip = false;
      return;
    }

    if (config.print) printf("gexon->right_gap = %d\n", gexon->right_gap);
    if (gexon->right_gap > 0) {
      tid_data.has_right_clip = false;
      return;
    }
    
    // Start position includes both insertion and soft clip
    int start_pos = seq_len - n_right_ins - n_right_clip;
    uint32_t total_bases = n_right_ins + n_right_clip;

    // Get pointer to the encoded sequence
    // uint8_t *seq = bam_get_seq(b);

    // Extract the bases (insertion + soft clip)
    char rescue_seq[total_bases + 1];
    for (int i = 0; i < total_bases; i++) {
      rescue_seq[i] = seq_nt16_str[bam_seqi(seq, start_pos + i)];
    }
    rescue_seq[total_bases] = '\0';  // Null terminate

    std::string remaining_qseq = rescue_seq;
    std::shared_ptr<GuideExon> current_gexon = gexon;
    
    const int MAX_EXONS = 20;
    int exons_processed = 0;

    int total_query_consumed = 0;
    while (!remaining_qseq.empty() && exons_processed < MAX_EXONS) {
      exons_processed++;
      
      // Get next guide exon
      std::shared_ptr<GuideExon> right_gexon;
      if (strand == '+') {
        if (!current_gexon->has_next) {
          if (exons_processed == 1) {
            if (config.print) printf("No exon to the right\n");
            tid_data.has_right_clip = false;
            return;
          } else {
            break;
          }
        }
        right_gexon = g2t->getGuideExonForTid(refid, strand, tid, 
          current_gexon->next_start, current_gexon->next_end);
      } else {
        if (!current_gexon->has_prev) {
          if (exons_processed == 1) {
            if (config.print) printf("No exon to the right\n");
            tid_data.has_right_clip = false;
            return;
          } else {
            break;
          }
        }
        right_gexon = g2t->getGuideExonForTid(refid, strand, tid, 
          current_gexon->prev_start, current_gexon->prev_end);
      }

      std::string next_gseq = right_gexon->seq;

      if (config.print) {
        auto tid_string = g2t->getTidName(tid);
        printf("\nGetting right alignment for tid = %s (exon %d)\n", 
               tid_string.c_str(), exons_processed);
        printf("remaining_qseq = %s (len=%zu)\n", remaining_qseq.c_str(), remaining_qseq.length());
        printf("next_gseq = %s (len=%zu)\n", next_gseq.c_str(), next_gseq.length());
        printf("curr_gexon_start = %d, curr_gexon_end = %d\n", current_gexon->start, current_gexon->end);
        printf("next_gexon_start = %d, next_gexon_end = %d\n", right_gexon->start, right_gexon->end);
      }

      AlignmentResult result = smith_waterman(remaining_qseq, next_gseq, START);
        // returns INCLUSIVE end values

      if (config.print) {
        printf("\nRIGHT CLIP ALIGNMENT RESULT:\n");
        printf("score = %d\n", result.score);
        printf("result.accuracy = %f\n", result.accuracy);
        printf("query_start = %d, query_end = %d\n", result.start_i, result.end_i);
        printf("ref_start = %d, ref_end = %d\n", result.start_j, result.end_j);
      }

      if (exons_processed == 1 && (result.accuracy <= 0.70  // 0.50 and 0 for 94 pearson
        || result.score < 10)) {
        tid_data.has_right_clip = false;
        if (config.print) printf("Score too low, giving up\n");
        return;
      }

      // Query exon coordinates based on where alignment occurs in guide
      uint32_t qstart = right_gexon->start + result.start_j;
      uint32_t qend = right_gexon->start + result.end_j;
      
      right_gexon->pos = result.pos + right_gexon->pos_start;

      if (exons_processed == 1) {
        // If we have an insertion, correct the current segment
        if (n_right_ins > 0) {
          segment.qexon.end = gexon->end; // remove insertion from this exon
        }
      }

      Segment right_clip;
      right_clip.is_valid = true;
      right_clip.qexon = GSeg(qstart, qend);
      right_clip.gexon = right_gexon;
      right_clip.has_qexon = true;
      right_clip.has_gexon = true;
      right_clip.status = RIGHTC_EXON;
      right_clip.is_small_exon = (remaining_qseq.length() <= config.small_exon_size);
      right_clip.cigar = result.cigar;
      right_clip.score = result.score;
      tid_data.segments.emplace_back(right_clip);
      tid_data.has_right_clip = true;

      // Calculate remaining query sequence
      int query_consumed = result.end_i + 1;

      // Always update remaining_qseq first
      if (query_consumed < remaining_qseq.length()) {
        remaining_qseq = remaining_qseq.substr(query_consumed);
      } else {
        remaining_qseq.clear();  // Everything consumed
      }

      // Then check if we should continue
      if (remaining_qseq.empty()) {
        break;
      }

      // Update for next iteration
      current_gexon = right_gexon;
      break;
    }

    // After while loop ends
    if (!remaining_qseq.empty()) {
      // Add final soft clip for any leftover bases
      tid_data.segments.back().cigar.add_operation(remaining_qseq.length(), BAM_CLIP_OVERRIDE);
    }
  }

  void ReadEvaluator::create_match(std::unordered_map<tid_t, TidData> &data, 
                                  std::shared_ptr<GuideExon> gexon, tid_t tid,
                                  char strand, char read_strand, bool has_gexon, 
                                  uint32_t ins_exon_len, ReadEvaluationConfig config) {
    ExonChainMatch match;
    match.align.fwpos = gexon->pos;
    match.align.rcpos = gexon->pos;
    match.transcript_len = gexon->transcript_len;
    match.align.strand = strand;
    match.align.cigar = std::make_shared<Cigar>();
    match.align.is_reverse = (strand != read_strand);
    match.align.similarity_score = 0;
    match.total_coverage = 0;
    match.total_operations = ins_exon_len;
    match.ref_consumed = 0;
    match.prev_op = BAM_CMATCH;
    match.junc_hits = 0;

    data[tid].match = match;
  }

  void ReadEvaluator::build_cigar_match(Segment &segment, TidData &tid_data,
                                        ExonChainMatch &match, bool first_match,
                                        bool last_match, ReadEvaluationConfig &config) {
    uint32_t qstart = segment.qexon.start;
    uint32_t qend = segment.qexon.end;
    uint32_t gstart = segment.gexon->start;
    uint32_t gend = segment.gexon->end;

    std::shared_ptr<Cigar> cigar = match.align.cigar;

    // Handle start boundary
    if (qstart < gstart) {
      uint32_t overhang = gstart - qstart;
      if (segment.status == FIRST_EXON || segment.status == ONLY_EXON) {
        if (!tid_data.has_left_clip) {
          cigar->add_operation(overhang, BAM_CSOFT_CLIP);
          match.total_operations += overhang;
          match.prev_op = BAM_CSOFT_CLIP;
          if (config.print) printf("Start of exon: adding softclip len %d\n", overhang);
        } else {
          //match.total_operations += overhang; // add penalty for this
        }
        
      } else if (segment.status == MIDDLE_EXON || segment.status == LAST_EXON
        || tid_data.has_left_clip) {
        cigar->add_operation(overhang, BAM_CINS);
        match.total_operations += overhang;
        if (match.prev_op == BAM_CDEL) {
          match.total_coverage += overhang;
        } else if (match.prev_op == BAM_CINS) {
          match.total_operations += (match.total_operations * 0.2); // penalty for two I in a row
        }
        match.prev_op = BAM_CINS;
        if (config.print) printf("Start of exon: adding insertion len %d\n", overhang);
      }
    } else if (gstart < qstart) {
      if (!first_match && (segment.status == MIDDLE_EXON || segment.status == LAST_EXON
        || tid_data.has_left_clip)) {
        uint32_t overhang = qstart - gstart;
        cigar->add_operation(overhang, BAM_CDEL);
        match.total_operations += overhang;
        match.ref_consumed += overhang;
        if (match.prev_op == BAM_CINS) {
          match.total_coverage += overhang;
        } else if (match.prev_op == BAM_CDEL) {
          match.total_operations += (match.total_operations * 0.2); // penalty for two D in a row
        }
        match.prev_op = BAM_CDEL;
        if (config.print) printf("Start of exon: adding deletion len %d\n", overhang);
      }
    // } else if (!tid_data.has_left_clip && first_match) {
    //   match.junc_hits++;
    // }
    } else {
      match.junc_hits++;
    }

    // Add match for overlapping region
    uint32_t overlap_start = std::max(qstart, gstart);
    uint32_t overlap_end = std::min(qend, gend);  
    if (overlap_end >= overlap_start) {
      uint32_t match_length = overlap_end - overlap_start;
      cigar->add_operation(match_length, BAM_CMATCH);
      match.total_operations += match_length;
      match.total_coverage += match_length;
      match.ref_consumed += match_length;
      match.prev_op = BAM_CMATCH;
      if (config.print) printf("Adding match len %d\n", match_length);
    }

    // Handle end boundary
    if (gend < qend) {
      uint32_t overhang = qend - gend;
      if (segment.status == LAST_EXON || segment.status == ONLY_EXON) {
        if (!tid_data.has_right_clip) {
          cigar->add_operation(overhang, BAM_CSOFT_CLIP);
          match.total_operations += overhang;
          match.prev_op = BAM_CSOFT_CLIP;
          if (config.print) printf("End of exon: adding softclip len %d\n", overhang);
        } else {
         // match.total_operations += overhang; // add penalty for this
        }
        
      } else if (segment.status == FIRST_EXON || segment.status == MIDDLE_EXON
        || tid_data.has_right_clip) {
        cigar->add_operation(overhang, BAM_CINS); 
        match.total_operations += overhang;
        if (match.prev_op == BAM_CDEL) {
          match.total_coverage += overhang;
        }
        match.prev_op = BAM_CINS;
        if (config.print) printf("End of exon: adding insertion len %d\n", overhang);
      }
    } else if (qend < gend) {
      if (!last_match && (segment.status == FIRST_EXON || segment.status == MIDDLE_EXON
        || tid_data.has_right_clip)) {
        uint32_t overhang = gend - qend;
        cigar->add_operation(overhang, BAM_CDEL); 
        match.total_operations += overhang;
        match.ref_consumed += overhang;
        if (match.prev_op == BAM_CINS) {
          match.total_coverage += overhang;
        }
        match.prev_op = BAM_CDEL;
        if (config.print) printf("End of exon: adding deletion len %d\n", overhang);
      }
    // } else if (!tid_data.has_right_clip && last_match){
    //   match.junc_hits++;
    // }
    } else {
      match.junc_hits++;
    }
  }

  void ReadEvaluator::build_cigar_ins(Segment &segment, TidData &tid_data, uint32_t k, uint32_t n,
                                      ExonChainMatch &match, ReadEvaluationConfig &config) {
    uint32_t qstart = segment.qexon.start;
    uint32_t qend = segment.qexon.end;
    uint32_t ins_length = qend - qstart;

    std::shared_ptr<Cigar> cigar = match.align.cigar;

    if (k == 0 || k == (n - 1)) {
      cigar->add_operation(ins_length, BAM_CSOFT_CLIP);
      match.prev_op = BAM_CSOFT_CLIP;
    } else {
      cigar->add_operation(ins_length, BAM_CINS);
      match.prev_op = BAM_CINS;
    }
    match.total_operations += ins_length;
    match.total_coverage += ins_length;
    if (config.print) printf("adding ins len %d\n", ins_length);
  }

  void ReadEvaluator::build_cigar_gap(Segment &segment, TidData &tid_data,
                                      ExonChainMatch &match, ReadEvaluationConfig &config) {
    uint32_t gstart = segment.gexon->start;
    uint32_t gend = segment.gexon->end;
    uint32_t gap_length = gend - gstart;

    std::shared_ptr<Cigar> cigar = match.align.cigar;

    cigar->add_operation(gap_length, BAM_CDEL);
    match.prev_op = BAM_CDEL;
    match.total_operations += gap_length;
    match.total_coverage += gap_length;
    match.ref_consumed += gap_length;
    if (config.print) printf("adding gap len %d\n", gap_length);
  }

  void ReadEvaluator::build_cigar_clip(Segment &segment, TidData &tid_data,
                                      ExonChainMatch &match, ReadEvaluationConfig &config) {

    std::shared_ptr<Cigar> cigar = match.align.cigar;
    for (auto &cig : segment.cigar.cigar) {
      uint8_t op = cig & BAM_CIGAR_MASK;
      uint32_t len = cig >> BAM_CIGAR_SHIFT;
      cigar->add_operation(len, op);
      
      if (op == BAM_CMATCH_OVERRIDE 
        || op == BAM_CDEL_OVERRIDE) {
        match.ref_consumed += len;
      }
    }
    if (config.print) printf("Added clip cigar to match cigar\n");
    match.align.clip_score += segment.score;
  }

  int32_t get_rand(uint32_t x) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, x - 1);
    return dis(gen);
  }

  void 
  ReadEvaluator::filter_by_similarity(std::unordered_map<tid_t, ExonChainMatch> &matches, 
                                      std::shared_ptr<g2tTree> g2t, ReadEvaluationConfig config) {
    for (auto it = matches.begin(); it != matches.end(); ) {
      auto &match = it->second;

      double similarity = (match.total_operations > 0)
        ? (match.total_coverage / match.total_operations)
        : 0.0;

      std::string tid_string = g2t->getTidName(it->first);

      if (similarity > config.similarity_threshold) {
        double x = ((similarity - config.similarity_threshold) 
          / (1.0 - config.similarity_threshold));
        match.align.similarity_score = (x * x * static_cast<double>(match.junc_hits + 1)); // static_cast<double>(match.junc_hits + 1)
        
        if (config.print) {
          printf("x**3 = %f\n", x * x * x);
          printf("match.align.junc_hits + 1 = %f\n", static_cast<double>(match.junc_hits + 1));
          printf("tid = %s, similarity_score = %f\n", tid_string.c_str(), match.align.similarity_score);
        }

        ++it;
      } else {
        it = matches.erase(it);
      }
    }

    for (auto it = matches.begin(); it != matches.end(); ) {
      auto &match = it->second;
      std::string tid_string = g2t->getTidName(it->first);
      pos_t pos = match.align.strand == '+' ? 
        match.align.fwpos : match.align.rcpos;
      if (pos + match.ref_consumed > match.transcript_len) {
        printf("Error: read %s for transcript %s extends beyond transcript end\n",
          config.name.c_str(), tid_string.c_str());
        printf("transcript length = %d, reference consumed = %d, pos = %d\n",
          match.transcript_len, match.ref_consumed, pos);
        it = matches.erase(it);
      } else {
        ++it;
      }
    }

    if (matches.empty()) return;

    // Determine primary/secondary alignment status
    auto best_it = matches.end();
    double best_score = -std::numeric_limits<double>::infinity();
    int32_t hit_index = 1;
    int count_at_best = 0;

    for (auto it = matches.begin(); it != matches.end(); ++it) {
      auto& match = it->second;
      match.align.hit_index = hit_index++;

      if (match.align.similarity_score > best_score) {
        best_score = match.align.similarity_score;
        best_it = it;
        count_at_best = 1;
      } else if (match.align.similarity_score == best_score) {
        count_at_best++;
      }
    }

    if (best_it != matches.end() && count_at_best == 1) {
      best_it->second.align.primary_alignment = true;
    } else {
      int32_t rand_idx = get_rand(matches.size());
      auto it = matches.begin();
      std::advance(it, rand_idx);
      it->second.align.primary_alignment = true;
    }
  }

  void print_segment(const Segment &s, int idx = -1) {
    if (idx >= 0) printf("  [%d] ", idx);
    printf("SEG valid=%d status=%d small=%d | ",
          s.is_valid, s.status, s.is_small_exon);

    if (s.has_qexon) {
      printf("Q[%d-%d] ",
            s.qexon.start, s.qexon.end);
    } else {
      printf("Q[.] ");
    }

    if (s.has_gexon && s.gexon) {
      printf("G[id=%d %d-%d]",
            s.gexon->exon_id,
            s.gexon->start,
            s.gexon->end);
    } else {
      printf("G[.]");
    }

    printf("\n");
  }

  std::unordered_map<tid_t, ExonChainMatch> 
  ReadEvaluator::evaluate_exon_chains(CReadAln *read, 
                                      read_id_t id, std::shared_ptr<g2tTree> g2t, 
                                      ReadEvaluationConfig config, uint8_t *seq,
                                      int seq_len) {
    uint32_t exon_count = read->segs.Count();
    refid_t refid = read->refid;
    config.name = read->brec->name();

    std::unordered_map<tid_t, ExonChainMatch> matches_by_strand;

    uint32_t n_left_ins = 0;
    uint32_t n_right_ins = 0;

    bool has_left_clip = false;
    bool has_right_clip = false;
    uint32_t n_left_clip = 0;
    uint32_t n_right_clip = 0;

    bool failure = false;

    if (config.name == "NM_021109_166_aligned_3877359_F_1_377_4") {
      config.print = true;
    } else {
      config.print = false;
    }
    //printf("read name = %s\n", config.name.c_str());
    //config.print = true;
    //GMessage("read name = %s\n", config.name.c_str());

    if (LONG_READS) {
      get_clips(read, config, failure, has_left_clip, has_right_clip,
        n_left_clip, n_right_clip);
    }

    std::unordered_map<tid_t, TidData> data;

    auto strands_to_check = get_strands_to_check(read);
    for (char strand : strands_to_check) {

      data.clear();

      Segment start_ignore;
      Segment end_ignore;

      failure = false;
      bool has_variants = false;
      bool has_gap = false;
      for (uint32_t j = 0; j < exon_count; j++) {
        get_intervals(data, read, j, exon_count, config, g2t,
          refid, strand, has_left_clip, has_right_clip, 
          start_ignore, end_ignore, has_gap, failure);

        //if (!config.print) failure = true;
        if (failure) break;

      } // end of exon loop

      // try next strand
      if (failure) continue;

      if (config.print) {
        printf("\n=== BEFORE CORRECTION ===\n");
        for (auto &pair : data) {
          tid_t tid = pair.first;
          auto &tid_data = pair.second;

          if (tid_data.elim) continue;

          printf("TID %d (%s), elim=%d\n",
                tid,
                g2t->getTidName(tid).c_str(),
                tid_data.elim);

          int k = 0;
          for (auto &s : tid_data.segments) {
            print_segment(s, k++);
          }
        }
        printf("\n");
      }

      // Correct for missed guide exons
      for (auto &pair : data) {
        tid_t tid = pair.first;
        auto &tid_data = pair.second;

        if (tid_data.elim) continue;

        correct_for_gaps(tid_data, tid, config, g2t,
          strand, refid);
      }

      if (config.print) {
        printf("\n=== AFTER CORRECTION ===\n");
        for (auto &pair : data) {
          tid_t tid = pair.first;
          auto &tid_data = pair.second;

          if (tid_data.elim) continue;

          printf("TID %d (%s), elim=%d\n",
                tid,
                g2t->getTidName(tid).c_str(),
                tid_data.elim);

          int k = 0;
          for (auto &s : tid_data.segments) {
            print_segment(s, k++);
          }
        }
        printf("\n");
      }

      // Add soft clip segments for each tid
      if (LONG_READS) {
        for (auto &pair : data) {
          tid_t tid = pair.first;
          auto &tid_data = pair.second;

          if (tid_data.elim) continue;

          // if (start_ignore.is_valid) {
          //   if (config.print) printf("start ignore is valid\n");
          //   //tid_data.segments.emplace(tid_data.segments.begin(), start_ignore);
          //   n_left_clip += (start_ignore.qexon.end - start_ignore.qexon.start);
          // }

          // if (end_ignore.is_valid) {
          //   if (config.print) printf("end ignore is valid\n");
          //   //tid_data.segments.emplace_back(end_ignore);
          //   n_right_clip += (end_ignore.qexon.end - end_ignore.qexon.start);
          // }

          if (tid_data.has_left_ins && USE_FASTA) {
            if (config.print) printf("has_left_ins\n");
            if (n_left_clip >= 5) {
              left_clip_rescue(tid_data, start_ignore, strand, g2t, refid, 
                tid, n_left_clip, tid_data.n_left_ins, config, read, 
                seq, seq_len);
            } else {
              tid_data.has_left_clip = false;
            }
          }
          else if (tid_data.has_left_clip && USE_FASTA) {
            if (n_left_clip >= 5) {
              left_clip_rescue(tid_data, start_ignore, strand, g2t, refid,
                tid, n_left_clip, 0, config, read, seq, seq_len);
            } else {
              tid_data.has_left_clip = false;
            }
          }

          if (tid_data.has_right_ins && USE_FASTA) {
            if (config.print) printf("has_right_ins\n");
            if (n_right_clip >= 5) {
              right_clip_rescue(tid_data, end_ignore, strand, g2t, refid,
                tid, n_right_clip, tid_data.n_right_ins, config, read, 
                seq, seq_len);
            } else {
              tid_data.has_right_clip = false;
            }
          }
          else if (tid_data.has_right_clip && USE_FASTA) {
            if (n_right_clip >= 5) {
              right_clip_rescue(tid_data, end_ignore, strand, g2t, refid,
                tid, n_right_clip, 0, config, read, seq, seq_len);
            } else {
              tid_data.has_right_clip = false;
            }
          }
          
        }

        if (config.print) {
          printf("\n=== AFTER CLIP RESCUE ===\n");
          for (auto &pair : data) {
            tid_t tid = pair.first;
            auto &tid_data = pair.second;

            if (tid_data.elim) continue;

            printf("TID %d (%s), elim=%d\n",
                  tid,
                  g2t->getTidName(tid).c_str(),
                  tid_data.elim);

            int k = 0;
            for (auto &s : tid_data.segments) {
              print_segment(s, k++);
            }
          }
          printf("\n");
        }
      }

      // Calculate matches using segments
      // segments should now reflect ignored query exons, ignored guide exons,
      // and rectified soft clips on either side
      for (auto &pair : data) {
        tid_t tid = pair.first;
        auto &tid_data = pair.second;

        if (tid_data.elim) continue;

        auto n_segments = tid_data.segments.size();
        uint32_t total_exon_matches = exon_count;

        std::string tid_string = g2t->getTidName(tid);

        if (config.print) printf("\nGet Matches, new tid = %d, %s\n", tid, tid_string.c_str());

        //int k_match = 0;
        uint32_t n_ins_exon = 0;
        bool match_created = false;

        uint32_t prev_s = 0;
        uint32_t prev_e = 0;

        int first_match_idx = -1;
        int last_match_idx = -1;

        for (auto &segment: tid_data.segments) {

          // if first exon is an INS_EXON, match starts next
          // keep track of ins_exon length to count operations correctly
          if (!match_created && segment.status == INS_EXON) {
            n_ins_exon += (segment.qexon.end - segment.qexon.start);
            first_match_idx++;
            last_match_idx++;
            continue;
          }

          // dunno what this is doing
          if (segment.has_gexon) {
            if (segment.gexon->start == prev_s 
              && segment.gexon->end == prev_e) {
                tid_data.elim = true;
                break;
              }
          }

          // Create exon chain matches 
          if (!match_created && segment.status != INS_EXON) {
            create_match(data, segment.gexon, tid, strand, read->strand, 
              segment.has_gexon, n_ins_exon, config);
            match_created = true;
            first_match_idx++;
            last_match_idx++;
            if (config.print) printf("match created, fwpos = %d, rcpos = %d\n", tid_data.match.align.fwpos, tid_data.match.align.rcpos);
          }

          else if (match_created // && segment.status != GAP_EXON
            && segment.status != INS_EXON) {
            // update rcpos
            last_match_idx++;
            if (strand == '-') {
              tid_data.match.align.rcpos = segment.gexon->pos;
              if (config.print) printf("Update rcpos: tid = %s, rcpos = %d\n", tid_string.c_str(), tid_data.match.align.rcpos);
            }
          }
        }

        auto& match = tid_data.match;
        if (!match.align.cigar) {
          if (config.print) printf("match cigar is null, tid = %d\n", tid);
          tid_data.elim = true;
        }

        int k_build = -1;
        bool first_match = false;
        bool last_match = false;

        for (auto &segment: tid_data.segments) {
          k_build++;

          first_match = (k_build == first_match_idx);
          last_match = (k_build == last_match_idx);

          if (tid_data.elim) continue;
          
          if (config.print) {
            print_segment(segment, k_build);
            printf("has_left_clip = %d, has_right_clip = %d\n", tid_data.has_left_clip, tid_data.has_right_clip);
            printf("is first match ? %d is last match ? %d\n", first_match, last_match);
            printf("first match idx %d last match idx %d\n", first_match_idx, last_match_idx);
          }

          // Build CIGAR as we process each exon
          bool is_match = (segment.status == FIRST_EXON || segment.status == MIDDLE_EXON
            || segment.status == LAST_EXON || segment.status == ONLY_EXON);
          bool is_ins = (segment.status == INS_EXON);
          bool is_gap = (segment.status == GAP_EXON);
          bool is_clip = (segment.status == LEFTC_EXON || segment.status == RIGHTC_EXON);

          if (is_match) {
            build_cigar_match(segment, tid_data, match, first_match, 
              last_match, config);
            
          } else if (is_ins) {
            build_cigar_ins(segment, tid_data, k_build, n_segments,
              match, config);
            if (k_build == 0 || k_build == (n_segments - 1)) {
              match.junc_hits -= 1;
            } else {
              match.junc_hits -= 2;
            }

          } else if (is_gap) {
            build_cigar_gap(segment, tid_data, match, config);
            match.junc_hits -= 2;
            first_match = false;

          } else if (is_clip) {
            build_cigar_clip(segment, tid_data, match, config);
            first_match = false;
          }
          
        } // end of segment loop

        // shouldn't happen
        if (match.junc_hits < 0) {
          match.junc_hits = 0;
        }

        // Record match for this strand
        if (!tid_data.elim) matches_by_strand[tid] = match;
      }

      if (config.print) {
        for (auto &pair : matches_by_strand) {
          std::string tid_string = g2t->getTidName(pair.first);
          printf("tid = %d, %s\n", pair.first, tid_string.c_str());
          printf("fwpos = %d\n", pair.second.align.fwpos);
          printf("rcpos = %d\n", pair.second.align.rcpos);
          printf("IDEAL CIGAR: ");
          for (const auto &cig : pair.second.align.cigar->cigar) {
            printf("%u%c ", cig >> BAM_CIGAR_SHIFT, 
              "MIDNSHP=XB,./;"[cig & BAM_CIGAR_MASK]);
          }
          printf("\n\n");
        }
      }

    } // end of strands to check

    if (!matches_by_strand.empty()) {
      filter_by_similarity(matches_by_strand, g2t, config);
    } else {
      if (failure)
        dropped_reads++;
      else unresolved_reads++;
    }
    return matches_by_strand;
  }

  std::unordered_map<tid_t, ExonChainMatch> 
  ShortReadEvaluator::evaluate(CReadAln * read, 
                              read_id_t id, 
                              std::shared_ptr<g2tTree> g2t,
                              uint8_t *seq,
                              int seq_len) {

    uint32_t max_clip = 25;
    uint32_t max_ins = 5;
    uint32_t max_gap = 5;
    uint32_t max_junc_gap = 0;
    float similarity_threshold = 0.80;

    ReadEvaluationConfig config = {
      max_clip,                       // max clip size
      max_ins,                        // max insertion to intervals
        // for when there exist no guides to explain a portion of an exon
      max_gap,                        // max gap in reference to intervals
      false,                          // ignore small exons?
      0,                              // small exon size
      max_junc_gap,                   // max junction gap
      similarity_threshold,
      false                           // print debug statements
    }; 

    auto matches = evaluate_exon_chains(read, id, g2t, config, seq, seq_len);
    return matches;
  }

  std::unordered_map<tid_t, ExonChainMatch> 
  LongReadEvaluator::evaluate(CReadAln * read, 
                              read_id_t id, 
                              std::shared_ptr<g2tTree> g2t,
                              uint8_t *seq,
                              int seq_len) {

    uint32_t max_clip = 40;
    uint32_t max_ins = 40;
    uint32_t max_gap = 40;
    uint32_t max_junc_gap = 40;
    float similarity_threshold = 0.60;
    
    ReadEvaluationConfig config = {
      max_clip,                       // max clip size
      max_ins,                        // max insertion to intervals
        // for when there exist no guides to explain a portion of an exon
      max_gap,                        // max gap in reference to intervals
      true,                           // ignore small exons?
      35,                             // small exon size
      max_junc_gap,                   // max junction gap
      similarity_threshold,
      false                           // print debug statements
    };   

    auto matches = evaluate_exon_chains(read, id, g2t, config, seq, seq_len); 
    return matches;
  }

} // namespace bramble