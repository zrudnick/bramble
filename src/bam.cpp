
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

namespace bramble {

  // Expand CIGAR into single per-base operations
  void expand_cigars(uint32_t* real_cigar, uint32_t n_cigar,
                    const Cigar& ideal_cigar,
                    std::vector<char> &real_expanded,
                    std::vector<char> &ideal_expanded,
                    uint32_t real_front_soft_clip,
                    uint32_t real_front_hard_clip) {

    // ------ Real cigar

    for (uint32_t k = 0; k < n_cigar; k++) {
      uint32_t op = real_cigar[k] & BAM_CIGAR_MASK;
      uint32_t len = real_cigar[k] >> BAM_CIGAR_SHIFT;
      
      if (op == BAM_CREF_SKIP) continue; // skip introns
      
      char op_char = "MIDNSHP=XB"[op];
      for (uint32_t i = 0; i < len; i++) {
        real_expanded.push_back(op_char);
      }
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
    for (const auto& pair : ideal_cigar.cigar) {
      uint32_t len = pair.first;
      uint32_t op = pair.second;
      
      char op_char = "MIDNSHP=XB"[op];

      uint32_t i = 0;
      while (i < len) {
        // Pad with '_' where real has an 'I'
        if ((j + n_pad < real_expanded.size())
          && real_expanded[j+n_pad] == 'I') {
          ideal_expanded.push_back('_');
          j++;
        } else {
          ideal_expanded.push_back(op_char);
          i++; j++;
        }
      }
    }
  }

  // Rules for choosing ops
  char merge_ops(char real_op, char ideal_op) {
    if (real_op == 'I' && ideal_op == '_') {
      return '_';
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
                          uint32_t* new_n_cigar,
                          CigarMem& mem) {
    // todo: change to # unique ops in expanded
    uint32_t* new_cigar = mem.get_mem(expanded.size() + 2);
    uint32_t new_idx = 0;
    
    char current_op = expanded[0];
    uint32_t run_len = 1;
    
    auto add_op = [&](char op_char, uint32_t len) {
      if (len == 0) return;
      uint32_t op;
      switch(op_char) {
        case 'M': case '=': case 'X': op = BAM_CMATCH; break;
        case 'I': op = BAM_CINS; break;
        case 'D': op = BAM_CDEL; break;
        case 'S': op = BAM_CSOFT_CLIP; break;
        //case 'H': op = BAM_CHARD_CLIP; break;
        case 'H': return; // just do nothing ::cry::
        // todo: figure out why hard clips cause errors
        case '_': return;
        default: return;
      }
      new_cigar[new_idx++] = (len << BAM_CIGAR_SHIFT) | op;
    };
    
    for (size_t i = 1; i < expanded.size(); i++) {
      if (expanded[i] == '_') {
        // do nothing
      }
      else if (expanded[i] == current_op) {
        run_len++;
      } else {
        add_op(current_op, run_len);
        current_op = expanded[i];
        run_len = 1;
      }
    }
    add_op(current_op, run_len);
    
    *new_n_cigar = new_idx;
    return new_cigar;
  }

  uint32_t* get_new_cigar(uint32_t* real_cigar, uint32_t n_real_cigar,
                        const Cigar& ideal_cigar, uint32_t* new_n_cigar, 
                        CigarMem& mem) {
    
    // Get front clip length from real cigar
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
    uint32_t* result = compress_cigar(merged, new_n_cigar, mem);

    // // Print debug info
    //   fprintf(stderr, "\nCIGAR SEQUENCE LENGTH CHANGED: %u -> %u\n", real_total_len, new_total_len);
    //   fprintf(stderr, "REAL CIGAR: ");
    //   for (uint32_t k = 0; k < n_real_cigar; k++) {
    //       uint32_t op = real_cigar[k] & BAM_CIGAR_MASK;
    //       uint32_t len = real_cigar[k] >> BAM_CIGAR_SHIFT;
    //       fprintf(stderr, "%u%c", len, "MIDNSHP=XB"[op]);
    //   }
    //   fprintf(stderr, "\nIDEAL CIGAR: ");
    //   for (const auto& pair : ideal_cigar.cigar) {
    //       fprintf(stderr, "%u%c ", pair.first, "MIDNSHP=XB"[pair.second]);
    //   }
    //   fprintf(stderr, "\nREAL EXPANDED:  ");
    //   for (char c : real_expanded) fprintf(stderr, "%c", c);
    //   fprintf(stderr, "\nIDEAL EXPANDED: ");
    //   for (char c : ideal_expanded) fprintf(stderr, "%c", c);
    //   fprintf(stderr, "\nMERGED:         ");
    //   for (char c : merged) fprintf(stderr, "%c", c);
    //   fprintf(stderr, "\nNEW CIGAR: ");
    //   for (uint32_t k = 0; k < *new_n_cigar; k++) {
    //       uint32_t op = result[k] & BAM_CIGAR_MASK;
    //       uint32_t len = result[k] >> BAM_CIGAR_SHIFT;
    //       fprintf(stderr, "%u%c", len, "MIDNSHP=XB"[op]);
    //   }
    //   fprintf(stderr, "\n--------------------\n");
    
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
                    CigarMem& mem, const Cigar& this_cigar) {

    uint32_t* new_cigar = nullptr;
    uint32_t new_n_cigar = 0;

    new_cigar = get_new_cigar(cigar, n_cigar, this_cigar, 
      &new_n_cigar, mem);

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
      pos = (int32_t)this_pair->r_pos;
    } else {
      tid = (int32_t)this_pair->m_tid;
      pos = (int32_t)this_pair->m_pos;
    }
    
    // Set basic paired flags
    b->core.flag |= BAM_FPAIRED;
    
    if (first_read && this_pair->read2->is_reverse) {
      b->core.flag |= BAM_FMREVERSE;
    } else if (!first_read && this_pair->read1->is_reverse) {
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

  void set_xs_tag(bam1_t* b, char xs_a) {
    uint8_t* xs = bam_aux_get(b, "XS");
    if (xs) bam_aux_del(b, xs);
    uint8_t new_xs = (uint8_t) xs_a;
    bam_aux_append(b, "XS", 'A', 1, &new_xs);
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