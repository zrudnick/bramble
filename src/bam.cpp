
#include <string>
#include <vector>
#include <algorithm>
#include <cstdlib>

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

  char merge_ops(char real_op, char ideal_op) {

    if ((real_op == BAM_CMATCH || real_op == BAM_CSOFT_CLIP) 
      && ideal_op == BAM_CLIP_OVERRIDE) {
      return BAM_CSOFT_CLIP;
    }
    if ((real_op == BAM_CMATCH || real_op == BAM_CSOFT_CLIP) 
      && ideal_op == BAM_CMATCH_OVERRIDE) {
      return BAM_CMATCH;
    }
    if ((real_op == BAM_CMATCH || real_op == BAM_CSOFT_CLIP) 
      && ideal_op == BAM_CINS_OVERRIDE) {
      return BAM_CINS;
    }
    if ((real_op == BAM_CMATCH || real_op == BAM_CSOFT_CLIP) 
      && ideal_op == BAM_CDEL_OVERRIDE) {
      return BAM_CDEL;
    }

    if (real_op == BAM_CDEL && ideal_op == BAM_CMATCH_OVERRIDE) {
      return BAM_CDEL;
    }

    if (real_op == BAM_CINS && ideal_op == BAM_CLIP_OVERRIDE) {
      return BAM_CSOFT_CLIP;
    }
    if (real_op == BAM_CINS && ideal_op == BAM_CMATCH_OVERRIDE) {
      return BAM_CINS;
    }

    if (ideal_op == BAM_CLIP_OVERRIDE) {
      return BAM_CSOFT_CLIP;
    }
    if (ideal_op == BAM_CMATCH_OVERRIDE) {
      return BAM_CMATCH;
    }
    if (ideal_op == BAM_CINS_OVERRIDE) {
      return BAM_CINS;
    }
    if (ideal_op == BAM_CDEL_OVERRIDE) { 
      // this takes care of the '_' and '.' case, rest may be unnecessary
      return BAM_CDEL;
    }

    if (real_op == BAM_CPAD) {
      return ideal_op;
    }
    
    if (real_op == BAM_CHARD_CLIP) {
      return BAM_CHARD_CLIP;
    }
    
    // Real I and ideal S - use S
    if (real_op == BAM_CINS && ideal_op == BAM_CSOFT_CLIP) {
      return BAM_CSOFT_CLIP;
    }
    
    // If ideal has S, D, or I, use it
    if (ideal_op == BAM_CSOFT_CLIP || ideal_op == BAM_CDEL 
      || ideal_op == BAM_CINS) {
      return ideal_op;
    }
    
    // If real has S, D, or I, use it
    if (real_op == BAM_CSOFT_CLIP || real_op == BAM_CDEL 
      || real_op == BAM_CINS) {
      return real_op;
    }
    
    // Matches (M, =, X all become M)
    if (ideal_op == BAM_CMATCH || ideal_op == BAM_CEQUAL) {
      return BAM_CMATCH;
    }
    if (ideal_op == BAM_CDIFF) {
      return BAM_CDIFF;
    }
    if (real_op == BAM_CMATCH || real_op == BAM_CEQUAL) {
      return BAM_CMATCH;
    }
    if (real_op == BAM_CDIFF) {
      return BAM_CDIFF;
    }
    
    return ideal_op; // fallback
  }

  uint32_t* merge_cigars(uint32_t* real_cigar, uint32_t n_real_cigar,
                        const Cigar& ideal_cigar, uint32_t* new_n_cigar,
                        uint32_t real_front_hard_clip,
                        uint32_t real_front_soft_clip, CigarMem &mem) {
    
    uint32_t n_ideal_cigar = ideal_cigar.n_cigar;
    uint32_t max_size = n_real_cigar + n_ideal_cigar;
    uint32_t* result = mem.get_mem(max_size);
    uint32_t result_idx = 0;
    
    uint32_t ri = 0;
    uint32_t ii = 0;
    uint32_t real_pos = 0;
    uint32_t ideal_pos = 0;
    
    auto add_op = [&](uint8_t op, uint32_t len) {
      if (len == 0) return;
      if (result_idx > 0 && bam_cigar_op(result[result_idx - 1]) == op) {
        result[result_idx - 1] += (len << 4);
      } else {
        result[result_idx++] = bam_cigar_gen(len, op);
      }
    };
    
    auto get_remaining = [](uint32_t cigar_val, uint32_t pos) {
      return bam_cigar_oplen(cigar_val) - pos;
    };
    
    // Front hard clips
    uint32_t clips_remaining = real_front_hard_clip;
    while (clips_remaining > 0 && ri < n_real_cigar) {
      uint32_t available = get_remaining(real_cigar[ri], real_pos);
      uint32_t chunk = (clips_remaining < available) ? clips_remaining : available;
      
      add_op(bam_cigar_op(real_cigar[ri]), chunk);
      
      clips_remaining -= chunk;
      real_pos += chunk;
      if (real_pos >= bam_cigar_oplen(real_cigar[ri])) {
        ri++;
        real_pos = 0;
      }
    }

    // Front soft clips
    clips_remaining = real_front_soft_clip;
    while (clips_remaining > 0 && ri < n_real_cigar) {
      uint8_t real_op = bam_cigar_op(real_cigar[ri]);
      uint8_t ideal_op = (ii < n_ideal_cigar) ? 
        bam_cigar_op(ideal_cigar.cigar[ii]) : 0xff;
      
      uint32_t real_remaining = get_remaining(real_cigar[ri], real_pos);
      uint32_t ideal_remaining = (ii < n_ideal_cigar) ? 
        get_remaining(ideal_cigar.cigar[ii], ideal_pos) : UINT32_MAX;
      
      bool is_override = (ii < n_ideal_cigar && 
        (ideal_op == BAM_CMATCH_OVERRIDE || ideal_op == BAM_CDEL_OVERRIDE || 
        ideal_op == BAM_CINS_OVERRIDE || ideal_op == BAM_CLIP_OVERRIDE));
      
      if (is_override) {
        if (ideal_op == BAM_CDEL_OVERRIDE) {
          uint32_t chunk = ideal_remaining;
          
          add_op(merge_ops(real_op, ideal_op), chunk);
          
          ideal_pos += chunk;
          
          if (ideal_pos >= bam_cigar_oplen(ideal_cigar.cigar[ii])) {
            ii++;
            ideal_pos = 0;
          }
        } else {
          uint32_t chunk = clips_remaining;
          if (chunk > real_remaining) chunk = real_remaining;
          if (chunk > ideal_remaining) chunk = ideal_remaining;
          
          add_op(merge_ops(real_op, ideal_op), chunk);
          
          clips_remaining -= chunk;
          real_pos += chunk;
          ideal_pos += chunk;
          
          if (real_pos >= bam_cigar_oplen(real_cigar[ri])) {
            ri++;
            real_pos = 0;
          }
          if (ideal_pos >= bam_cigar_oplen(ideal_cigar.cigar[ii])) {
            ii++;
            ideal_pos = 0;
          }
        }
      } else {
        uint32_t chunk = clips_remaining;
        if (chunk > real_remaining) chunk = real_remaining;
        
        add_op(merge_ops(real_op, ideal_op), chunk);
        
        clips_remaining -= chunk;
        real_pos += chunk;
        
        if (real_pos >= bam_cigar_oplen(real_cigar[ri])) {
          ri++;
          real_pos = 0;
        }
      }
    }
    
    // Main loop
    while (ri < n_real_cigar || ii < n_ideal_cigar) {
      if (ri >= n_real_cigar) {
        uint32_t remaining = get_remaining(ideal_cigar.cigar[ii], ideal_pos);
        add_op(bam_cigar_op(ideal_cigar.cigar[ii]), remaining);
        ii++;
        ideal_pos = 0;
        continue;
      }
        
      if (ii >= n_ideal_cigar) {
        uint32_t remaining = get_remaining(real_cigar[ri], real_pos);
        add_op(bam_cigar_op(real_cigar[ri]), remaining);
        ri++;
        real_pos = 0;
        continue;
      }
        
      uint8_t real_op = bam_cigar_op(real_cigar[ri]);
      uint8_t ideal_op = bam_cigar_op(ideal_cigar.cigar[ii]);
      uint32_t real_remaining = get_remaining(real_cigar[ri], real_pos);
      uint32_t ideal_remaining = get_remaining(ideal_cigar.cigar[ii], ideal_pos);
      
      if (real_op == BAM_CREF_SKIP) {
        ri++;
        real_pos = 0;
      } else if (real_op == BAM_CDEL && 
        (ideal_op == BAM_CSOFT_CLIP || ideal_op == BAM_CLIP_OVERRIDE ||
        ideal_op == BAM_CINS || ideal_op == BAM_CINS_OVERRIDE)) {
        uint32_t chunk = (real_remaining < ideal_remaining) ? 
          real_remaining : ideal_remaining;
        
        real_pos += chunk;
        ideal_pos += chunk;
        
        if (real_pos >= bam_cigar_oplen(real_cigar[ri])) {
          ri++;
          real_pos = 0;
        }
        if (ideal_pos >= bam_cigar_oplen(ideal_cigar.cigar[ii])) {
          ii++;
          ideal_pos = 0;
        }
      } else if (real_op == BAM_CINS) {
        add_op(BAM_CINS, real_remaining);
        ri++;
        real_pos = 0;
      } else if (ideal_op == BAM_CDEL || ideal_op == BAM_CDEL_OVERRIDE) {
        add_op(BAM_CDEL, ideal_remaining);
        ii++;
        ideal_pos = 0;
      } else {
        uint32_t chunk = (real_remaining < ideal_remaining) ? 
          real_remaining : ideal_remaining;
        uint8_t merged_op = merge_ops(real_op, ideal_op);
        add_op(merged_op, chunk);
        
        real_pos += chunk;
        ideal_pos += chunk;
        
        if (real_pos >= bam_cigar_oplen(real_cigar[ri])) {
          ri++;
          real_pos = 0;
        }
        if (ideal_pos >= bam_cigar_oplen(ideal_cigar.cigar[ii])) {
          ii++;
          ideal_pos = 0;
        }
      }
    }
    
    *new_n_cigar = result_idx;
    return result;
  }

  void print_debug(uint32_t* real_cigar, uint32_t n_real_cigar,
                  const Cigar& ideal_cigar, uint32_t* new_n_cigar,
                  uint32_t* result) {
    fprintf(stderr, "REAL CIGAR: ");
    for (uint32_t k = 0; k < n_real_cigar; k++) {
      uint32_t op = real_cigar[k] & BAM_CIGAR_MASK;
      uint32_t len = real_cigar[k] >> BAM_CIGAR_SHIFT;
      fprintf(stderr, "%u%c", len, "MIDNSHP=XB,./;"[op]);
    }

    fprintf(stderr, "\nIDEAL CIGAR: ");
    for (int i = 0; i < ideal_cigar.n_cigar; ++i) {
      auto cig = ideal_cigar.cigar[i];
      fprintf(stderr, "%u%c ", cig >> BAM_CIGAR_SHIFT, 
        "MIDNSHP=XB,./;"[cig & BAM_CIGAR_MASK]);
    }

    std::string real_expanded, ideal_expanded, merged_expanded;
    size_t ri = 0, ii = 0;
    uint32_t real_pos = 0, ideal_pos = 0;

    while (ri < n_real_cigar || ii < ideal_cigar.n_cigar) {
      uint8_t real_op = 0, ideal_op = 0;
      
      if (ri < n_real_cigar) {
        real_op = real_cigar[ri] & BAM_CIGAR_MASK;
      }
      if (ii < ideal_cigar.n_cigar) {
        ideal_op = ideal_cigar.cigar[ii] & BAM_CIGAR_MASK;
      }
      
      if (ri >= n_real_cigar) {
        real_expanded += '_';
        ideal_expanded += "MIDNSHP=XB,./;"[ideal_op];
        merged_expanded += "MIDNSHP=XB,./;"[ideal_op];
        ideal_pos++;
        if (ideal_pos >= (ideal_cigar.cigar[ii] >> BAM_CIGAR_SHIFT)) {
          ii++;
          ideal_pos = 0;
        }
      } else if (ii >= ideal_cigar.n_cigar) {
        real_expanded += "MIDNSHP=XB,./;"[real_op];
        ideal_expanded += '_';
        merged_expanded += "MIDNSHP=XB,./;"[real_op];
        real_pos++;
        if (real_pos >= (real_cigar[ri] >> BAM_CIGAR_SHIFT)) {
          ri++;
          real_pos = 0;
        }
      } else if (real_op == BAM_CREF_SKIP) {
        real_pos++;
        if (real_pos >= (real_cigar[ri] >> BAM_CIGAR_SHIFT)) {
          ri++;
          real_pos = 0;
        }
      } else if (real_op == BAM_CHARD_CLIP) {
        real_expanded += 'H';
        ideal_expanded += '_';
        merged_expanded += 'H';
        real_pos++;
        if (real_pos >= (real_cigar[ri] >> BAM_CIGAR_SHIFT)) {
          ri++;
          real_pos = 0;
        }
      } else if (real_op == BAM_CSOFT_CLIP) {
        if (ideal_op == BAM_CMATCH_OVERRIDE) {
          real_expanded += 'S';
          ideal_expanded += ',';
          merged_expanded += 'M';
          ideal_pos++;
          real_pos++;
        } else if (ideal_op == BAM_CDEL_OVERRIDE) {
          real_expanded += '_';
          ideal_expanded += '.';
          merged_expanded += 'D';
          ideal_pos++;
        } else if (ideal_op == BAM_CINS_OVERRIDE) {
          real_expanded += 'S';
          ideal_expanded += '/';
          merged_expanded += 'I';
          ideal_pos++;
          real_pos++;
        } else if (ideal_op == BAM_CLIP_OVERRIDE) {
          real_expanded += 'S';
          ideal_expanded += ';';
          merged_expanded += 'S';
          ideal_pos++;
          real_pos++;
        } else {
          real_expanded += 'S';
          ideal_expanded += '_';
          merged_expanded += 'S';
          real_pos++;
        }
        if (real_pos >= (real_cigar[ri] >> BAM_CIGAR_SHIFT)) {
          ri++;
          real_pos = 0;
        }
        if (ideal_pos >= (ideal_cigar.cigar[ii] >> BAM_CIGAR_SHIFT)) {
          ii++;
          ideal_pos = 0;
        }
      } else if (real_op == BAM_CINS) {
        real_expanded += 'I';
        ideal_expanded += '_';
        merged_expanded += 'I';
        real_pos++;
        if (real_pos >= (real_cigar[ri] >> BAM_CIGAR_SHIFT)) {
          ri++;
          real_pos = 0;
        }
      } else if (ideal_op == BAM_CDEL || ideal_op == BAM_CDEL_OVERRIDE) {
        real_expanded += '_';
        ideal_expanded += "MIDNSHP=XB,./;"[ideal_op];
        merged_expanded += 'D';
        ideal_pos++;
        if (ideal_pos >= (ideal_cigar.cigar[ii] >> BAM_CIGAR_SHIFT)) {
          ii++;
          ideal_pos = 0;
        }
      } else {
        real_expanded += "MIDNSHP=XB,./;"[real_op];
        ideal_expanded += "MIDNSHP=XB,./;"[ideal_op];
        merged_expanded += "MIDNSHP=XB,./;"[merge_ops(real_op, ideal_op)];
        real_pos++;
        ideal_pos++;
        if (real_pos >= (real_cigar[ri] >> BAM_CIGAR_SHIFT)) {
          ri++;
          real_pos = 0;
        }
        if (ideal_pos >= (ideal_cigar.cigar[ii] >> BAM_CIGAR_SHIFT)) {
          ii++;
          ideal_pos = 0;
        }
      }
    }

    fprintf(stderr, "\nREAL EXPANDED:  %s", real_expanded.c_str());
    fprintf(stderr, "\nIDEAL EXPANDED: %s", ideal_expanded.c_str());
    fprintf(stderr, "\nMERGED:         %s", merged_expanded.c_str());

    fprintf(stderr, "\nNEW CIGAR: ");
    for (uint32_t k = 0; k < *new_n_cigar; k++) {
      uint32_t op = result[k] & BAM_CIGAR_MASK;
      uint32_t len = result[k] >> BAM_CIGAR_SHIFT;
      fprintf(stderr, "%u%c", len, "MIDNSHP=XB,./;"[op]);
    }
    fprintf(stderr, "\n--------------------\n");
  }

  uint32_t* get_new_cigar(uint32_t* real_cigar, uint32_t n_real_cigar,
                        const Cigar& ideal_cigar, uint32_t* new_n_cigar,
                        CigarMem& mem, int32_t& nm) {

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

    uint32_t* result = merge_cigars(real_cigar, n_real_cigar,
      ideal_cigar, new_n_cigar, real_front_hard_clip, 
      real_front_soft_clip, mem);

    // debug cigar strings
    //print_debug(real_cigar, n_real_cigar, ideal_cigar, new_n_cigar, result);
    
    return result;
  }

  // Copy CIGAR memory for CIGAR updates
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
      int32_t r_pos = (this_pair->r_align.strand == '+') ? this_pair->r_align.fwpos : this_pair->r_align.rcpos;
      int32_t m_pos = (this_pair->m_align.strand == '+') ? this_pair->m_align.fwpos : this_pair->m_align.rcpos;

      int32_t my_pos   = first_read ? r_pos : m_pos;
      int32_t mate_pos = first_read ? m_pos : r_pos;


      b->core.mtid = b->core.tid;
      b->core.mpos = mate_pos;
      b->core.flag |= BAM_FPROPER_PAIR;
      
      // SAM spec: positive for leftmost read, negative for rightmost
      if (my_pos <= mate_pos) {
        b->core.isize = (mate_pos + b->core.l_qseq) - my_pos;
      } else {
        b->core.isize = -((my_pos + b->core.l_qseq) - mate_pos);
      }
              
    // Mates map to different transcripts
    } else {
      int32_t r_pos = (this_pair->r_align.strand == '+') ? this_pair->r_align.fwpos : this_pair->r_align.rcpos;
      int32_t m_pos = (this_pair->m_align.strand == '+') ? this_pair->m_align.fwpos : this_pair->m_align.rcpos;
      
      b->core.mtid = first_read ? this_pair->m_tid : this_pair->r_tid;
      b->core.mpos = first_read ? m_pos : r_pos;  // use the same fwpos/rcpos logic
      b->core.isize = 0;
      b->core.flag &= ~BAM_FPROPER_PAIR;
    }
  }

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

  void set_ts_tag(bam1_t* b, char ts_a) {
    uint8_t* ts = bam_aux_get(b, "ts");
    if (ts) bam_aux_del(b, ts);
    // uint8_t new_ts = (uint8_t) ts_a;
    // bam_aux_append(b, "ts", 'A', 1, &new_ts);
  }

  void set_as_tag(bam1_t* b, AlignInfo &align, uint8_t qual) {
    uint8_t* as = bam_aux_get(b, "AS");
    int32_t gn_as = 0;
    if (as) gn_as = bam_aux2i(as);
    if (as) bam_aux_del(b, as);

    double score;
    double x = align.similarity_score;
    double y = static_cast<double>(align.clip_score);
    score = (static_cast<double>(gn_as) + y) * x;
    int32_t new_as = static_cast<int32_t>(score);
    bam_aux_append(b, "AS", 'i', sizeof(int32_t), (uint8_t*)&new_as);
  }

  int reverse_complement_bam(bam1_t *b) {
    // validate record
    if (!b) return -1;
    if (b->core.l_qseq <= 0) {
      // No sequence to reverse-complement, but still need to reverse CIGAR
      // and flip the flag so the alignment is correct in transcript coordinates.
      uint32_t n_cigar = b->core.n_cigar;
      uint32_t *cigar = bam_get_cigar(b);
      for (uint32_t i = 0; i < n_cigar / 2; i++) {
        uint32_t temp = cigar[i];
        cigar[i] = cigar[n_cigar - 1 - i];
        cigar[n_cigar - 1 - i] = temp;
      }
      b->core.flag ^= BAM_FREVERSE;
      return 0;
    }

    int len = b->core.l_qseq;

    uint8_t *seq = bam_get_seq(b);
    uint8_t *qual = bam_get_qual(b);

    static const uint8_t comp_table[16] = {
      15, // 0 (unused)
      8,  // 1 A -> T
      4,  // 2 C -> G
      15, // 3
      2,  // 4 G -> C
      15,15,15,
      1,  // 8 T -> A
      15,15,15,15,15,15,15
    };

    // Allocate buffer for the reversed complemented sequence
    // Each byte holds two bases
    std::vector<uint8_t> tmp((len + 1) / 2);
    for (int i = 0; i < len; ++i) {
      uint8_t nt = bam_seqi(seq, len - 1 - i);
      uint8_t nt_complement = comp_table[nt];
      bam_set_seqi(tmp.data(), i, nt_complement);
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