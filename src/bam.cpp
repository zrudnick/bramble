
#include <set>
#include <unordered_set>
#include <string>
#include <vector>
#include <algorithm>
#include <memory>
#include <cstdlib>

#include "types.h"
#include "bundles.h"
#include "bramble.h"
#include "g2t.h"
#include "evaluate.h"
#include "bam.h"
#include "htslib/sam.h"

#ifndef NOTHREADS
#include "GThreads.h"
#endif

int seen_s;
extern double similarity_threshold;
extern bool LONG_READS;

namespace bramble {

  // Expand CIGAR into single per-base operations
  void expand_cigars(uint32_t* real_cigar, uint32_t n_cigar,
                    const Cigar& ideal_cigar,
                    std::vector<char> &real_expanded,
                    std::vector<char> &ideal_expanded,
                    uint32_t real_front_soft_clip,
                    uint32_t real_front_hard_clip) {

    // ------ Real cigar

    uint req = 0;
    for (uint32_t k = 0; k < n_cigar; k++) {
      uint32_t op = real_cigar[k] & BAM_CIGAR_MASK;
      uint32_t len = real_cigar[k] >> BAM_CIGAR_SHIFT;
      
      if (op == BAM_CREF_SKIP) continue; // skip introns
      
      char op_char = "MIDNSHP=XB"[op];
      for (uint32_t i = 0; i < len; i++) {
        real_expanded.push_back(op_char);
      }
      if (k == 0 && op_char == 'H') req += 1;
      if (k > req && op_char == 'S') seen_s += len;
      if (k > req && op_char == 'H') seen_s += len;
    }

    // ------ Ideal cigar

    // Pad ideal with '_' to match real's front hard clip
    for (uint32_t i = 0; i < real_front_hard_clip; i++) {
      ideal_expanded.push_back('_');
    }

    // Pad ideal with '_' to match real's front soft clip
    for (uint32_t i = 0; i < real_front_soft_clip; i++) {
      ideal_expanded.push_back('_');
    }

    uint32_t n_pad = real_front_hard_clip + real_front_soft_clip;
    uint32_t j = 0;
    uint32_t n_d = 0;
    for (const auto& pair : ideal_cigar.cigar) {
      uint32_t len = pair.first;
      uint32_t op = pair.second;
      
      char op_char = "MIDNSHP=XB"[op];

      uint32_t i = 0;
      while (i < len) {
        // Pad with '_' where real has an 'I'
        if (((j+n_pad-n_d) < real_expanded.size())
          && real_expanded[j+n_pad-n_d] == 'I') {
          ideal_expanded.push_back('_');
          j++;
        } else if (op_char == 'D') {
          ideal_expanded.push_back(op_char);
          i++; j++; n_d++;
        } else {
          ideal_expanded.push_back(op_char);
          i++; j++;
        }
      }
    }

    // Pad real_expanded with '_' where ideal has 'D'
    if (n_d > 0) {
      uint32_t j = 0;
      for (uint32_t i = 0; i < ideal_expanded.size(); i++) {
        if (ideal_expanded[i] == 'D') {
          real_expanded.insert(real_expanded.begin() + j, '_');
        }
        j++;
      }
    }
  }

  // Rules for choosing ops
  char merge_ops(char real_op, char ideal_op) {
    if (real_op == 'I' && ideal_op == '_') {
      return 'I';
    }

    // Added for end soft/hard clips
    if (ideal_op == '*') {
      return real_op;
    }
    // Unused
    if (real_op == '*') {
      return ideal_op;
    }
    if (real_op == 'H') {
      return 'H';
    }
    // If real D and ideal S, 
    // add '_' to avoid adding extra sequence
    // allow everything on the left to shift right
    // because whatever, it doens't matter anymore
    if (real_op == 'D' && ideal_op == 'S') {
      return '_';
    }
    if (real_op == 'I' && ideal_op == 'S') {
      return 'S';
    }
    if (real_op == 'D' && ideal_op == 'I') {
      return '_';
    }
    // If ideal has S, D, or I, use it
    if (ideal_op == 'S' || ideal_op == 'D' || ideal_op == 'I') {
      return ideal_op;
    }
    // If real has S, D, or I, use it
    if (real_op == 'S' || real_op == 'D' || real_op == 'I') {
      return real_op;
    }
    // Otherwise use M (or =/X)
    if (ideal_op == 'M' || ideal_op == '=' || ideal_op == 'X') {
      return 'M';
    }
    if (real_op == 'M' || real_op == '=' || real_op == 'X') {
      return 'M';
    }
    // Fallback
    return real_op != '*' ? real_op : ideal_op;
  }

  // Compress expanded CIGAR back to standard format
  uint32_t* compress_cigar(const std::vector<char>& expanded, 
                          uint32_t* new_n_cigar, CigarMem& mem, int32_t &nm) {
    if (expanded.empty()) {
      *new_n_cigar = 0;
      return nullptr;
    }
    
    std::vector<char> processed = expanded;
    size_t size = processed.size();
    uint32_t num_runs = 0;
    char prev_op = '\0';
    
    for (size_t i = 0; i < size; ) {
      char curr = processed[i];
      
      if (curr == '_') {
        i++;
        continue;
      }
      
      if (curr == 'S') {
        size_t s_start = i;
        while (i < size && processed[i] == 'S') i++;
        
        if (i < size && processed[i] == 'I') {
          size_t i_start = i;
          while (i < size && processed[i] == 'I') i++;
          
          if (i < size && processed[i] == 'S') {
            for (size_t j = i_start; j < i; j++) {
              processed[j] = 'S';
            }
          }
        }
        
        if (prev_op != 'S') {
          num_runs++;
          prev_op = 'S';
        }
      } else {
        if (curr != prev_op) {
          num_runs++;
          prev_op = curr;
        }
        i++;
      }
    }
    
    uint32_t* new_cigar = mem.get_mem(expanded.size());
    uint32_t new_idx = 0;
    
    char current_op = processed[0];
    uint32_t run_len = 1;
    
    auto add_op = [&](char op_char, uint32_t len) {
      if (len == 0) return;
      uint32_t op;
      switch(op_char) {
        case 'M': case '=': case 'X': op = BAM_CMATCH; break;
        case 'I': op = BAM_CINS; nm++; break;
        case 'D': op = BAM_CDEL; nm++; break;
        case 'S': op = BAM_CSOFT_CLIP; break;
        case 'H': return;
        case '_': return;
        default: nm++; return;
      }
      new_cigar[new_idx++] = (len << BAM_CIGAR_SHIFT) | op;
    };
    
    for (size_t i = 1; i < size; i++) {
      if (processed[i] == '_') {
        // do nothing
      }
      else if (processed[i] == current_op) {
        run_len++;
      } else {
        add_op(current_op, run_len);
        current_op = processed[i];
        run_len = 1;
      }
    }
    add_op(current_op, run_len);
    
    *new_n_cigar = new_idx;
    return new_cigar;
  }

  uint32_t* get_new_cigar(uint32_t* real_cigar, uint32_t n_real_cigar,
                        const Cigar& ideal_cigar, uint32_t* new_n_cigar, 
                        CigarMem& mem, int32_t &nm) {
    
    // Get front clip length from real cigar
    uint32_t real_front_hard_clip = 0;
    uint32_t real_front_soft_clip = 0;
    uint32_t cigar_idx = 0;
    seen_s = 0;

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

    // Expand both cigars
    std::vector<char> real_expanded;
    std::vector<char> ideal_expanded;
    expand_cigars(real_cigar, n_real_cigar, ideal_cigar,
      real_expanded, ideal_expanded, real_front_soft_clip,
      real_front_hard_clip);
    
    // Merge position by position
    size_t max_len = std::max(real_expanded.size(), ideal_expanded.size());
    std::vector<char> merged;
    
    for (size_t i = 0; i < max_len; i++) {
      char real_op = (i < real_expanded.size()) ? real_expanded[i] : '*';
      char ideal_op = (i < ideal_expanded.size()) ? ideal_expanded[i] : '*';
      merged.push_back(merge_ops(real_op, ideal_op));
    }
    
    // Compress back to CIGAR
    uint32_t* result = compress_cigar(merged, new_n_cigar, mem, nm);

    bool ideal_has_s = false;
    for (const auto& pair : ideal_cigar.cigar) {
      uint32_t len = pair.first;
      uint32_t op = pair.second;
      if (op == BAM_CSOFT_CLIP && len > 100) ideal_has_s = true;
    }

    bool wrong_len = (real_expanded.size() != ideal_expanded.size() + seen_s);

    // Print debug info
    // if (ideal_has_s) {
      // fprintf(stderr, "REAL CIGAR: ");
      // for (uint32_t k = 0; k < n_real_cigar; k++) {
      //     uint32_t op = real_cigar[k] & BAM_CIGAR_MASK;
      //     uint32_t len = real_cigar[k] >> BAM_CIGAR_SHIFT;
      //     fprintf(stderr, "%u%c", len, "MIDNSHP=XB"[op]);
      // }
      // fprintf(stderr, "\nIDEAL CIGAR: ");
      // for (const auto& pair : ideal_cigar.cigar) {
      //     fprintf(stderr, "%u%c ", pair.first, "MIDNSHP=XB"[pair.second]);
      // }
      // fprintf(stderr, "\nREAL EXPANDED:  ");
      // for (char c : real_expanded) fprintf(stderr, "%c", c);
      // fprintf(stderr, "\nIDEAL EXPANDED: ");
      // for (char c : ideal_expanded) fprintf(stderr, "%c", c);
      // fprintf(stderr, "\nMERGED:         ");
      // for (char c : merged) fprintf(stderr, "%c", c);
      // fprintf(stderr, "\nNEW CIGAR: ");
      // for (uint32_t k = 0; k < *new_n_cigar; k++) {
      //     uint32_t op = result[k] & BAM_CIGAR_MASK;
      //     uint32_t len = result[k] >> BAM_CIGAR_SHIFT;
      //     fprintf(stderr, "%u%c", len, "MIDNSHP=XB"[op]);
      // }
      // fprintf(stderr, "\n--------------------\n");
    // }
    
    return result;
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
    int32_t tid;
    int32_t pos;
    if (first_read) {
      tid = (int32_t)this_pair->r_tid;
      pos = (int32_t)this_pair->r_align.pos;
    } else {
      tid = (int32_t)this_pair->m_tid;
      pos = (int32_t)this_pair->m_align.pos;
    }
    
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
      // todo: replace read_size with correct measure
      int32_t isize = 0;
      if (first_read) {
        isize = ((b->core.mpos + this_pair->read2->read_size) - b->core.pos);
      } else {
        isize = -((b->core.pos + this_pair->read2->read_size) - b->core.mpos); // negative for second read in pair
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
    // uint8_t new_xs = (uint8_t) xs_a;
    // bam_aux_append(b, "XS", 'A', 1, &new_xs);
  }

  void set_as_tag(bam1_t* b, double similarity_score) {
    uint8_t* as = bam_aux_get(b, "AS");
    if (as) bam_aux_del(b, as);
    double score = std::pow(1.0+(similarity_score - similarity_threshold), 3.0)*100.0;
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

} // namespace bramble