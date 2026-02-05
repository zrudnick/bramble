
#include <set>
#include <unordered_set>
#include <string>
#include <vector>
#include <algorithm>
#include <memory>
#include <cstdlib>
#include <cmath>
#include <random>

#include "types.h"
#include "bramble.h"
#include "g2t.h"
#include "evaluate.h"
#include "bam.h"
#include "htslib/sam.h"

#ifndef NOTHREADS
#include "GThreads.h"
#endif

extern bool LONG_READS;

namespace bramble {

  // ============================================================================
  // STEP 1: Merge adjacent I and D into M
  // ============================================================================

  Cigar merge_indels(const Cigar& ideal_cigar) {
    Cigar result;
    if (ideal_cigar.cigar.empty()) return result;
    
    int32_t i_count = 0;
    int32_t d_count = 0;
    
    auto flush_indels = [&]() {
      if (i_count == 0 && d_count == 0) return;
      
      // Convert overlapping I and D to M
      int32_t overlap = std::min(i_count, d_count);
      if (overlap > 0) {
        result.add_operation((uint32_t)overlap, BAM_CMATCH);
        i_count -= overlap;
        d_count -= overlap;
      }
      
      // Add remaining I or D
      if (i_count > 0) {
        result.add_operation((uint32_t)i_count, BAM_CINS);
        i_count = 0;
      }
      if (d_count > 0) {
        result.add_operation((uint32_t)d_count, BAM_CDEL);
        d_count = 0;
      }
    };
    
    for (const auto &cig : ideal_cigar.cigar) {
      uint8_t op = cig & BAM_CIGAR_MASK;
      uint32_t len = cig >> BAM_CIGAR_SHIFT;
      
      if (op == BAM_CINS) {
        i_count += len;
      } else if (op == BAM_CDEL) {
        d_count += len;
      } else {
        flush_indels();
        result.add_operation(len, op);
      }
    }
    
    flush_indels();
    return result;
  }

  // ============================================================================
  // STEP 2: Convert CIGAR to character arrays
  // ============================================================================

  // Expand packed BAM CIGAR (uint32_t array)
  void expand_cigar_simple(uint32_t* cigar, uint32_t n_cigar,
                          std::vector<char>& expanded) {
    for (uint32_t k = 0; k < n_cigar; k++) {
      uint8_t op = cigar[k] & BAM_CIGAR_MASK;
      uint32_t len = cigar[k] >> BAM_CIGAR_SHIFT;
      
      // Skip introns (N operations)
      if (op == BAM_CREF_SKIP) continue;
      
      char op_char = "MIDNSHP=XB,./;"[op];
      for (uint32_t i = 0; i < len; i++) {
        expanded.push_back(op_char);
      }
    }
  }

  // Expand Cigar struct (vector of pairs)
  void expand_cigar_struct(const Cigar& cigar, std::vector<char>& expanded) {
    for (const auto &cig : cigar.cigar) {
      uint8_t op = cig & BAM_CIGAR_MASK;
      uint32_t len = cig >> BAM_CIGAR_SHIFT;
      
      char op_char = "MIDNSHP=XB,./;"[op];
      for (uint32_t i = 0; i < len; i++) {
        expanded.push_back(op_char);
      }
    }
  }

  // ============================================================================
  // STEP 3: Pad ideal CIGAR to match real's leading clips
  // ============================================================================

  void pad_ideal_for_leading_clips(std::vector<char>& ideal_exp,
                                  uint32_t front_hard_clip,
                                  uint32_t front_soft_clip) {
    // Total padding needed at the front
    uint32_t total_padding = front_hard_clip + front_soft_clip;
    if (total_padding == 0) return;
    
    std::vector<char> padded;

    for (uint32_t i = 0; i < front_hard_clip; i++) {
      padded.push_back('_');
    }
    
    // Add padding at the front, but check if ideal has ,./ overrides
    // that should replace the padding at those positions
    for (uint32_t i = front_hard_clip; i < front_soft_clip; i++) {
      if (i < ideal_exp.size() && 
          (ideal_exp[i] == ',' || ideal_exp[i] == '.' 
          || ideal_exp[i] == '/' || ideal_exp[i] == ';')) {
        // do nothing
      } else {
        padded.push_back('_');
      }
    }
    
    for (uint32_t i = 0; i < ideal_exp.size(); i++) {
      padded.push_back(ideal_exp[i]);
    }
    
    ideal_exp = std::move(padded);
  }

  // ============================================================================
  // STEP 4: Align two expanded CIGARs by inserting padding characters
  // ============================================================================


  void align_expanded_cigars(std::vector<char>& real_exp,
                          std::vector<char>& ideal_exp) {
    // Build new aligned vectors instead of inserting (avoids iterator invalidation)
    std::vector<char> aligned_real;
    std::vector<char> aligned_ideal;
    
    size_t real_pos = 0;
    size_t ideal_pos = 0;
    
    // Safety counter to prevent infinite loops
    size_t max_iterations = (real_exp.size() + ideal_exp.size()) * 2;
    size_t iterations = 0;
    
    while (real_pos < real_exp.size() || ideal_pos < ideal_exp.size()) {
      // Safety check
      if (++iterations > max_iterations) {
        fprintf(stderr, "ERROR: align_expanded_cigars exceeded max iterations, possible infinite loop\n");
        break;
      }
      
      // Handle edge cases where one is exhausted
      if (real_pos >= real_exp.size()) {
        // Real is done, pad it to match remaining ideal
        aligned_real.push_back('_');
        aligned_ideal.push_back(ideal_exp[ideal_pos]);
        ideal_pos++;
        continue;
      }
      if (ideal_pos >= ideal_exp.size()) {
        // Ideal is done, pad it to match remaining real
        aligned_real.push_back(real_exp[real_pos]);
        aligned_ideal.push_back('_');
        real_pos++;
        continue;
      }
      
      char r = real_exp[real_pos];
      char i = ideal_exp[ideal_pos];
      
      // Special case: ideal has '.' (deletion override)
      // This means ideal wants to insert here, so pad real
      if (i == '.') {
        aligned_real.push_back('_');
        aligned_ideal.push_back(i);
        ideal_pos++;
        // Note: real_pos is NOT incremented
      }
      // Insertions in real = deletions in ideal's perspective
      // Add padding to ideal AND keep repeating the current ideal character
      else if (r == 'I') {
        aligned_real.push_back(r);
        aligned_ideal.push_back('_');
        real_pos++;
        // Note: ideal_pos is NOT incremented - same char gets used next iteration
      }
      // Deletions in ideal = insertions in real's perspective
      // Add padding to real AND keep repeating the current real character
      else if (i == 'D') {
        aligned_real.push_back('_');
        aligned_ideal.push_back(i);
        ideal_pos++;
        // Note: real_pos is NOT incremented - same char gets used next iteration
      }
      // Both advance normally
      else {
        aligned_real.push_back(r);
        aligned_ideal.push_back(i);
        real_pos++;
        ideal_pos++;
      }
    }
    
    // Replace original vectors with aligned versions
    real_exp = std::move(aligned_real);
    ideal_exp = std::move(aligned_ideal);
  }

  // ============================================================================
  // STEP 5: Merge operations with clear priority rules
  // ============================================================================

  char merge_ops(char real_op, char ideal_op) {
    // PRIORITY 1: Special case - real insertion with ideal padding
    if (real_op == 'I' && ideal_op == '_') {
      return 'I';
    }
    
    // PRIORITY 2: Special overrides - ideal's ,./ override real's soft clips
    // ; = soft clip override (should become S)
    // , = match override (should become M)
    // . = deletion override (should become D)
    // / = insertion override (should become I)
    if ((real_op == 'M' || real_op == 'S') && ideal_op == ';') {
      return 'S';
    }
    if ((real_op == 'M' || real_op == 'S') && ideal_op == ',') {
      return 'M';
    }
    if ((real_op == 'M' || real_op == 'S') && ideal_op == '/') {
      return 'I';
    }
    if ((real_op == 'M' || real_op == 'S') && ideal_op == '.') {
      return 'D';
    }

    if (real_op == 'D' && ideal_op == ';') {
      return '_';
    }
    if (real_op == 'D' && ideal_op == ',') {
      return 'D';
    }
    if (real_op == 'D' && ideal_op == '/') {
      return '_';
    }
    if (real_op == 'D' && ideal_op == '.') {
      return '_';   // doesn't happen
    }

    if (real_op == 'I' && ideal_op == ';') {
      return 'S';
    }
    if (real_op == 'I' && ideal_op == ',') {
      return 'I';
    }
    if (real_op == 'D' && ideal_op == '/') {
      return '_';
    }
    if (real_op == 'I' && ideal_op == '.') {
      return '_';   // doesn't happen
    }

    if (ideal_op == ';') {
      return 'S';
    }
    if (ideal_op == ',') {
      return 'M';
    }
    if (ideal_op == '/') {
      return 'I';
    }
    if (ideal_op == '.') { // this takes care of the '_' and '.' case, rest may be unnecessary
      return 'D';
    }
    
    // PRIORITY 3: Handle end soft/hard clips (padding markers)
    if (ideal_op == '*') {
      return real_op;
    }
    if (real_op == '*') {
      return ideal_op;
    }
    
    // PRIORITY 4: Hard clips always win
    if (real_op == 'H') {
      return 'H';
    }
    
    // PRIORITY 5: Special edge cases
    // Real D and ideal S - add padding to avoid extra sequence
    if (real_op == 'D' && ideal_op == 'S') {
      return '_';
    }
    // Real I and ideal S - use S
    if (real_op == 'I' && ideal_op == 'S') {
      return 'S';
    }
    // Real D and ideal I - use padding
    if (real_op == 'D' && ideal_op == 'I') {
      return '_';
    }
    
    // PRIORITY 6: If ideal has S, D, or I, use it
    if (ideal_op == 'S' || ideal_op == 'D' || ideal_op == 'I') {
      return ideal_op;
    }
    
    // PRIORITY 7: If real has S, D, or I, use it
    if (real_op == 'S' || real_op == 'D' || real_op == 'I') {
      return real_op;
    }
    
    // PRIORITY 8: Matches (M, =, X all become M)
    if (ideal_op == 'M' || ideal_op == '=' || ideal_op == 'X') {
      return 'M';
    }
    if (real_op == 'M' || real_op == '=' || real_op == 'X') {
      return 'M';
    }
    
    // PRIORITY 9: Padding from alignment
    if (real_op == '_') return ideal_op;
    if (ideal_op == '_') return real_op;
    
    // Fallback
    return real_op != '*' ? real_op : ideal_op;
  }

  // ============================================================================
  // STEP 5: Compress expanded CIGAR back to standard format
  // ============================================================================

  Cigar compress_cigar(const std::vector<char>& expanded, int32_t& nm) {
    Cigar result;
    
    if (expanded.empty()) {
      return result;
    }
    
    // First pass: clean up S-I-S patterns (convert middle I to S)
    std::vector<char> cleaned = expanded;
    for (size_t i = 1; i < cleaned.size() - 1; i++) {
      if (cleaned[i] == 'I') {
        // Look backward for S
        bool has_s_before = false;
        for (size_t j = i; j > 0; j--) {
          if (cleaned[j-1] == 'S') {
            has_s_before = true;
            break;
          }
          if (cleaned[j-1] != 'I') break;
        }
        
        // Look forward for S
        bool has_s_after = false;
        for (size_t j = i; j < cleaned.size() - 1; j++) {
          if (cleaned[j+1] == 'S') {
            has_s_after = true;
            break;
          }
          if (cleaned[j+1] != 'I') break;
        }
        
        // If I is sandwiched between S operations, convert to S
        if (has_s_before && has_s_after) {
          cleaned[i] = 'S';
        }
      }
    }
    
    // Run-length encode into Cigar struct
    size_t i = 0;
    while (i < cleaned.size()) {
      // Skip padding characters
      if (cleaned[i] == '_') {
        i++;
        continue;
      }
      
      char op_char = cleaned[i];
      uint32_t run_len = 0;
      
      // Count consecutive identical operations
      while (i < cleaned.size() && (cleaned[i] == op_char || cleaned[i] == '_')) {
        if (cleaned[i] == op_char) run_len++;
        i++;
      }
      
      if (run_len == 0) continue;
      
      // Convert to BAM operation code
      uint32_t op;
      switch(op_char) {
        case 'M': case '=': case 'X': op = BAM_CMATCH; break;
        case 'I': op = BAM_CINS; nm += run_len; break;
        case 'D': op = BAM_CDEL; nm += run_len; break;
        case 'S': op = BAM_CSOFT_CLIP; break;
        case 'H': op = BAM_CHARD_CLIP; break;
        case 'N': op = BAM_CREF_SKIP; break;
        case 'P': op = BAM_CPAD; break;
        default: continue; // Skip unknown operations
      }
      
      result.add_operation(run_len, op);
    }
    
    return result;
  }


  // ============================================================================
  // MAIN FUNCTION: Merge real and ideal CIGARs
  // ============================================================================

  uint32_t* get_new_cigar(uint32_t* real_cigar, uint32_t n_real_cigar,
                        const Cigar& ideal_cigar, uint32_t* new_n_cigar,
                        CigarMem& mem, int32_t& nm) {
    
    Cigar merged_ideal = merge_indels(ideal_cigar);

    // Expand both CIGARs to character arrays
    std::vector<char> real_expanded;
    std::vector<char> ideal_expanded;
    
    expand_cigar_simple(real_cigar, n_real_cigar, real_expanded);
    expand_cigar_struct(merged_ideal, ideal_expanded);  // Use struct version!

    uint32_t real_front_hard_clip = 0;
    uint32_t real_front_soft_clip = 0;
    uint32_t cigar_idx = 0;
    
    // Check for leading hard clip
    if (n_real_cigar > 0 && 
        ((real_cigar[0] & BAM_CIGAR_MASK) == BAM_CHARD_CLIP)) {
      real_front_hard_clip = real_cigar[0] >> BAM_CIGAR_SHIFT;
      cigar_idx++;
    }
    
    // Check for leading soft clip (after potential hard clip)
    if (cigar_idx < n_real_cigar && 
        ((real_cigar[cigar_idx] & BAM_CIGAR_MASK) == BAM_CSOFT_CLIP)) {
      real_front_soft_clip = real_cigar[cigar_idx] >> BAM_CIGAR_SHIFT;
    }

    pad_ideal_for_leading_clips(ideal_expanded, real_front_hard_clip, real_front_soft_clip);
    align_expanded_cigars(real_expanded, ideal_expanded);

    std::vector<char> merged;
    size_t max_len = std::max(real_expanded.size(), ideal_expanded.size());

    for (size_t i = 0; i < max_len; i++) {
      char real_op = (i < real_expanded.size()) ? real_expanded[i] : '*';
      char ideal_op = (i < ideal_expanded.size()) ? ideal_expanded[i] : '*';
      merged.push_back(merge_ops(real_op, ideal_op));
    }
    
    Cigar compressed = compress_cigar(merged, nm);
    Cigar final_cigar = merge_indels(compressed);
    
    uint32_t* final_result = mem.get_mem(final_cigar.cigar.size());
    for (size_t k = 0; k < final_cigar.cigar.size(); k++) {
      uint32_t len = final_cigar.cigar[k] >> BAM_CIGAR_SHIFT;
      uint8_t op = final_cigar.cigar[k] & BAM_CIGAR_MASK;
      final_result[k] = (len << BAM_CIGAR_SHIFT) | op;
    }
    *new_n_cigar = final_cigar.cigar.size();

    // Debug output
    // fprintf(stderr, "REAL CIGAR: ");
    // for (uint32_t k = 0; k < n_real_cigar; k++) {
    //   uint32_t op = real_cigar[k] & BAM_CIGAR_MASK;
    //   uint32_t len = real_cigar[k] >> BAM_CIGAR_SHIFT;
    //   fprintf(stderr, "%u%c", len, "MIDNSHP=XB,./;"[op]);
    // }
    // fprintf(stderr, "\nIDEAL CIGAR: ");
    // for (const auto& cig : ideal_cigar.cigar) {
    //   fprintf(stderr, "%u%c ", cig >> BAM_CIGAR_SHIFT, 
    //     "MIDNSHP=XB,./;"[cig & BAM_CIGAR_MASK]);
    // }
    // fprintf(stderr, "\nREAL EXPANDED:  ");
    // for (char c : real_expanded) fprintf(stderr, "%c", c);
    // fprintf(stderr, "\nIDEAL EXPANDED: ");
    // for (char c : ideal_expanded) fprintf(stderr, "%c", c);
    // fprintf(stderr, "\nMERGED:         ");
    // for (char c : merged) fprintf(stderr, "%c", c);
    // fprintf(stderr, "\nNEW CIGAR: ");
    // for (uint32_t k = 0; k < *new_n_cigar; k++) {
    //   uint32_t op = final_result[k] & BAM_CIGAR_MASK;
    //   uint32_t len = final_result[k] >> BAM_CIGAR_SHIFT;
    //   fprintf(stderr, "%u%c", len, "MIDNSHP=XB,./;"[op]);
    // }
    // fprintf(stderr, "\n--------------------\n");
    
    return final_result;
  }

  // Copy CIGAR memory for intron removal
  uint8_t* copy_cigar_memory(bam1_t* b, uint32_t new_n_cigar, 
                            uint32_t old_n_cigar, uint32_t new_m_data,
                            const uint32_t* new_cigar,
                            int l_qname, int l_qseq, int l_aux) {
    
    // Reallocate BAM data
    if (new_m_data > b->m_data) {
      b->m_data = new_m_data;
      b->data = (uint8_t*)realloc(b->data, b->m_data);
    }
    uint8_t* data = b->data;

    // BAM record layout = [QNAME][CIGAR][SEQ][QUAL][AUX]
    int old_seq_offset = l_qname + (old_n_cigar << 2);
    int new_seq_offset = l_qname + (new_n_cigar << 2);
    int seq_size = (l_qseq + 1) >> 1;
    int tail_size = seq_size + l_qseq + l_aux;

    memmove(data + new_seq_offset, data + old_seq_offset, tail_size);
    memcpy(data + l_qname, new_cigar, (size_t)new_n_cigar << 2);

    b->core.n_cigar = new_n_cigar;
    b->l_data = new_m_data;

    return data;
  }

  bool update_cigar(bam1_t* b, uint32_t* cigar, uint32_t n_cigar,
                    CigarMem& mem, const Cigar& ideal_cigar, int32_t &nm) {

    uint32_t* new_cigar = nullptr;
    uint32_t new_n_cigar = 0;

    new_cigar = get_new_cigar(cigar, n_cigar, ideal_cigar, 
      &new_n_cigar, mem, nm);

    // Calculate new data layout sizes
    // BAM record layout = [QNAME][CIGAR][SEQ][QUAL][AUX]
    int l_qname = b->core.l_qname;
    int l_qseq  = b->core.l_qseq;   // quality seq length
    int old_cigar_bytes = (n_cigar << 2);
    int new_cigar_bytes = (new_n_cigar << 2);
    int seq_size = (l_qseq + 1) >> 1;
    int l_aux = b->l_data - (l_qname + old_cigar_bytes + seq_size + l_qseq);

    // Calculate total data size (bytes) after replacing the CIGAR
    uint32_t new_m_data = (uint32_t)(l_qname + new_cigar_bytes + seq_size + l_qseq + l_aux);

    // Copy new cigar and shift seq/qual/aux
    copy_cigar_memory(b, new_n_cigar, n_cigar, new_m_data,
                      new_cigar, l_qname, l_qseq, l_aux);

    return true;
  }

  // Set mate information for a read
  void set_mate_info(bam1_t* b, BamInfo* this_pair, bool first_read) {

    // Handle unpaired reads
    if (!this_pair->is_paired) {
      // Clear all paired-end related flags
      b->core.flag &= ~(BAM_FPAIRED | BAM_FPROPER_PAIR | BAM_FMREVERSE);
      
      // Set mate reference to unmapped
      b->core.mtid = -1;  // -1 indicates unmapped mate
      b->core.mpos = -1;  // -1 indicates unmapped mate position
      b->core.isize = 0;  // No insert size for unpaired reads
      
      return;
    }
      
    // Get match tid and position
    int32_t tid = b->core.tid;
    int32_t pos = b->core.pos;
    
    // Set basic paired flags
    b->core.flag |= BAM_FPAIRED;
    
    bool read_is_reverse =  (this_pair->r_align.strand == '-');
    bool mate_is_reverse =  (this_pair->m_align.strand == '-');
    if (first_read && read_is_reverse) {
      b->core.flag |= BAM_FMREVERSE;
    } else if (!first_read && mate_is_reverse) {
      b->core.flag |= BAM_FMREVERSE;
    }
    
    // Mates map to same transcript
    if (this_pair->same_transcript) {
      // Both mates map to same transcript - use "=" for mate reference 
      b->core.mtid = tid; // Same as current read's tid
      b->core.mpos = pos;
      b->core.flag |= BAM_FPROPER_PAIR;
      
      // Calculate insert size for same transcript
      int32_t isize = 0;
      if (first_read) {
        if (b->core.mpos > b->core.pos) {
          isize = ((b->core.mpos + this_pair->read2->read_size) - b->core.pos);
        } else {
          isize = -((b->core.pos + this_pair->read1->read_size) - b->core.mpos);
        }
        
      } else {
        if (b->core.pos > b->core.mpos) {
          isize = -((b->core.pos + this_pair->read1->read_size) - b->core.mpos);
        } else {
          isize = ((b->core.mpos + this_pair->read2->read_size) - b->core.pos);
        }
        
      }
      b->core.isize = isize;
              
    // Mates map to different transcripts
    } else {
      b->core.mtid = tid;
      b->core.mpos = pos;
      b->core.isize = 0;
      b->core.flag |= BAM_FPROPER_PAIR;
    }
  }

  // Add new NH tag (after removing current one)
  void set_nh_tag(bam1_t* b, int32_t nh_i) {
    uint8_t* nh = bam_aux_get(b, "NH");
    if (nh) bam_aux_del(b, nh);
    bam_aux_append(b, "NH", 'i', 4, (uint8_t*)&nh_i);
  }

  void set_nm_tag(bam1_t* b, int32_t nm_i) {
    uint8_t* nm = bam_aux_get(b, "NM");
    if (nm) bam_aux_del(b, nm);
    bam_aux_append(b, "NM", 'i', 4, (uint8_t*)&nm_i);
  }

  void set_hi_tag(bam1_t* b, int32_t hi_i) {
    uint8_t* hi = bam_aux_get(b, "HI");
    if (hi) bam_aux_del(b, hi);
    bam_aux_append(b, "HI", 'i', 4, (uint8_t*)&hi_i);
  }

  void set_xs_tag(bam1_t* b, char xs_a) {
    uint8_t* xs = bam_aux_get(b, "XS");
    if (xs) bam_aux_del(b, xs);
    uint8_t new_xs = (uint8_t) xs_a;
    bam_aux_append(b, "XS", 'A', 1, &new_xs);
  }

  void set_as_tag(bam1_t* b, AlignInfo &align, uint8_t qual) {
    uint8_t* as = bam_aux_get(b, "AS");
    int32_t gn_as = 0;
    if (as) gn_as = bam_aux2i(as);
    if (as) bam_aux_del(b, as);

    double score;
    if (!LONG_READS) {
      score = std::pow(1.0 + align.similarity_score, 3.0) * 100.0;
    } else {
      double x = align.similarity_score;
      double y = static_cast<double>(align.clip_score);
      score = (static_cast<double>(gn_as) + y) * x;
    }
    int32_t new_as = static_cast<int32_t>(score);
    bam_aux_append(b, "AS", 'i', sizeof(int32_t), (uint8_t*)&new_as);
  }

  void remove_extra_tags(bam1_t* b) {
    uint8_t* sa = bam_aux_get(b, "SA");
    if (sa) bam_aux_del(b, sa);
    uint8_t* ms = bam_aux_get(b, "ms");
    if (ms) bam_aux_del(b, ms);
    uint8_t* nn = bam_aux_get(b, "nn");
    if (nn) bam_aux_del(b, nn);
    uint8_t* ts = bam_aux_get(b, "ts");
    if (ts) bam_aux_del(b, ts);
    uint8_t* tp = bam_aux_get(b, "tp");
    if (tp) bam_aux_del(b, tp);
    uint8_t* cm = bam_aux_get(b, "cm");
    if (cm) bam_aux_del(b, cm);
    uint8_t* s1 = bam_aux_get(b, "s1");
    if (s1) bam_aux_del(b, s1);
    uint8_t* s2 = bam_aux_get(b, "s2");
    if (s2) bam_aux_del(b, s2);
    uint8_t* de = bam_aux_get(b, "de");
    if (de) bam_aux_del(b, de);
    uint8_t* rl = bam_aux_get(b, "rl");
    if (rl) bam_aux_del(b, rl);
  }

  int reverse_complement_bam(bam1_t *b) {
    // validate record
    if (!b) return -1;
    if (b->core.l_qseq <= 0) return 0;

    int len = b->core.l_qseq;

    uint8_t *seq = bam_get_seq(b);
    uint8_t *qual = bam_get_qual(b);

    // Allocate buffer for the reversed complemented sequence
    // Each byte holds two bases
    std::vector<uint8_t> tmp((len + 1) / 2);
    for (int i = 0; i < len; ++i) {
      uint8_t nt = bam_seqi(seq, len - 1 - i);
      uint8_t nt_complement;
      switch (nt) {
        case 1: nt_complement = 8; break; // A -> T
        case 2: nt_complement = 4; break; // C -> G
        case 4: nt_complement = 2; break; // G -> C
        case 8: nt_complement = 1; break; // T -> A
        default: nt_complement = 15; break; // N or others -> N
      }
      bam_set_seqi(tmp, i, nt_complement);
    }

    memcpy(seq, tmp.data(), (len + 1) / 2);

    if (qual && qual[0] != 0xff) { // qualities exist - need to be reversed as well
      for (int i = 0; i < len / 2; ++i) {
        uint8_t qtmp = qual[i];
        qual[i] = qual[len - 1 - i];
        qual[len - 1 - i] = qtmp;
      }
    }

    // reverse cigar operations
    uint32_t n_cigar = b->core.n_cigar;
    uint32_t *cigar = bam_get_cigar(b);
    for (uint32_t i = 0; i < n_cigar / 2; i++) {
      uint32_t temp = cigar[i];
      cigar[i] = cigar[n_cigar - 1 - i];
      cigar[n_cigar - 1 - i] = temp;
    }

    // flip the reverse flag
    b->core.flag ^= BAM_FREVERSE;

    // TODO: need to handle MD (and potentially other tags) here, since invalid after reversal
    return 0;
  }

} // namespace bramble