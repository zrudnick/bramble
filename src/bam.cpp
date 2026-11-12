
#include <string>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <cctype>
#include <cstring>

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

    if (real_op == BAM_CDEL && (ideal_op == BAM_CSOFT_CLIP
      || ideal_op == BAM_CLIP_OVERRIDE)) {
      return '_';
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
      if (op == '_') return; // deletion is soft clipped
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
    
    for (uint32_t i = 1; i + 1 < result_idx; i++) {
      if (bam_cigar_op(result[i]) != BAM_CINS) continue;
      uint8_t prev = bam_cigar_op(result[i-1]);
      uint8_t next = bam_cigar_op(result[i+1]);
      if ((prev == BAM_CSOFT_CLIP || prev == BAM_CHARD_CLIP) &&
        (next == BAM_CSOFT_CLIP || next == BAM_CHARD_CLIP)) {
        result[i] = bam_cigar_gen(bam_cigar_oplen(result[i]), prev);
      }
    }

    uint32_t new_idx = 0;
    for (uint32_t i = 0; i < result_idx; i++) {
      uint8_t op = bam_cigar_op(result[i]);
      uint32_t len = bam_cigar_oplen(result[i]);
      if (new_idx > 0 && bam_cigar_op(result[new_idx-1]) == op)
        result[new_idx-1] += (len << 4);
      else
        result[new_idx++] = bam_cigar_gen(len, op);
    }
    result_idx = new_idx;

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
      } else if (ii >= ideal_cigar.n_cigar) {
        real_expanded += "MIDNSHP=XB,./;"[real_op];
        ideal_expanded += '_';
        merged_expanded += "MIDNSHP=XB,./;"[real_op];
        real_pos++;
      } else if (real_op == BAM_CREF_SKIP) {
        real_pos++;
      } else if (real_op == BAM_CHARD_CLIP) {
        real_expanded += 'H';
        ideal_expanded += '_';
        merged_expanded += 'H';
        real_pos++;
      } else if (real_op == BAM_CDEL && (ideal_op == BAM_CSOFT_CLIP
        || ideal_op == BAM_CLIP_OVERRIDE)) {
        real_expanded += 'D';
        ideal_expanded += "MIDNSHP=XB,./;"[ideal_op];
        merged_expanded += '_';
        ideal_pos++;
        real_pos++;
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
      } else if (real_op == BAM_CINS) {
        real_expanded += 'I';
        ideal_expanded += '_';
        merged_expanded += 'I';
        real_pos++;
      } else if (ideal_op == BAM_CDEL || ideal_op == BAM_CDEL_OVERRIDE) {
        real_expanded += '_';
        ideal_expanded += "MIDNSHP=XB,./;"[ideal_op];
        merged_expanded += 'D';
        ideal_pos++;
      } else {
        real_expanded += "MIDNSHP=XB,./;"[real_op];
        ideal_expanded += "MIDNSHP=XB,./;"[ideal_op];
        merged_expanded += "MIDNSHP=XB,./;"[merge_ops(real_op, ideal_op)];
        real_pos++;
        ideal_pos++;
      }

      if (real_pos >= (real_cigar[ri] >> BAM_CIGAR_SHIFT)) {
        ri++;
        real_pos = 0;
      }
      if (ideal_pos >= (ideal_cigar.cigar[ii] >> BAM_CIGAR_SHIFT)) {
        ii++;
        ideal_pos = 0;
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
    // print_debug(real_cigar, n_real_cigar, ideal_cigar, new_n_cigar, result);
    
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

  // Count the MD mismatches that fall in `keep` (sorted, disjoint query
  // intervals [start, end), genome orientation), walking the MD tag and the
  // genome CIGAR together; no allocations (this runs once per output record).
  // Returns -1 if MD is malformed or does not fit the CIGAR.
  static int64_t count_md_mismatches(const char* md, const std::vector<uint32_t>& gcigar,
                                     const std::vector<std::pair<uint32_t, uint32_t>>& keep) {
    size_t ci = 0;   // current genome CIGAR op
    uint32_t off = 0;  // aligned bases already used in op ci
    uint32_t q = 0;    // query position at the start of op ci
    auto is_aligned = [](int op) {
      return op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF;
    };
    // move to the next unused reference-aligned base; false if none left
    auto next_aligned = [&]() {
      while (ci < gcigar.size()) {
        int op = bam_cigar_op(gcigar[ci]);
        uint32_t len = bam_cigar_oplen(gcigar[ci]);
        if (is_aligned(op) && off < len) return true;
        if (bam_cigar_type(op) & 1) q += len;  // consumes query
        ci++;
        off = 0;
      }
      return false;
    };

    int64_t count = 0;
    size_t kp = 0;
    for (const char* p = md; *p; ) {
      if (isdigit((unsigned char)*p)) {
        uint64_t n = strtoul(p, (char**)&p, 10);
        while (n > 0) {  // matching bases
          if (!next_aligned()) return -1;
          uint32_t take = std::min<uint64_t>(n, bam_cigar_oplen(gcigar[ci]) - off);
          off += take;
          n -= take;
        }
      } else if (*p == '^') {  // deleted reference bases: no query position
        p++;
        while (*p && isalpha((unsigned char)*p)) p++;
      } else if (isalpha((unsigned char)*p)) {
        if (!next_aligned()) return -1;
        uint32_t pos = q + off++;
        while (kp < keep.size() && keep[kp].second <= pos) kp++;
        if (kp < keep.size() && keep[kp].first <= pos) count++;
        p++;
      } else {
        return -1;
      }
    }
    return next_aligned() ? -1 : count;  // MD must cover every aligned base
  }

  // NM for the projected alignment: mismatches of bases still aligned outside
  // clip rescue (from MD; else genome NM minus genome indels), plus mismatches
  // in rescued bases, plus every inserted/deleted base of the projected CIGAR.
  // MD itself is dropped later by keep_only_known_tags. Call after update_cigar
  // and before reverse_complement_bam; gcigar is the genome CIGAR.
  void set_projected_nm(bam1_t* b, const std::vector<uint32_t>& gcigar,
                        const AlignInfo& align) {
    uint8_t* nm_tag = bam_aux_get(b, "NM");
    uint8_t* md_tag = bam_aux_get(b, "MD");
    const char* md = md_tag ? bam_aux2Z(md_tag) : nullptr;
    if (!nm_tag && !md) return;  // nothing known about mismatches
    int64_t genome_nm = nm_tag ? bam_aux2i(nm_tag) : 0;

    // aligned (M/=/X) query intervals of the projected CIGAR, outside the
    // ranges realigned by clip rescue
    uint32_t* cigar = bam_get_cigar(b);
    uint32_t n_cigar = b->core.n_cigar;
    uint32_t qlen = bam_cigar2qlen(n_cigar, cigar);
    uint32_t lo = std::min(align.rescue_left, qlen);
    uint32_t hi = qlen - std::min(align.rescue_right, qlen);
    thread_local std::vector<std::pair<uint32_t, uint32_t>> keep;
    keep.clear();
    int64_t new_indels = 0;
    uint32_t q = 0;
    for (uint32_t i = 0; i < n_cigar; i++) {
      int op = bam_cigar_op(cigar[i]);
      uint32_t len = bam_cigar_oplen(cigar[i]);
      if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
        uint32_t s = std::max(q, lo), e = std::min(q + len, hi);
        if (s < e) keep.emplace_back(s, e);
      }
      if (op == BAM_CINS || op == BAM_CDEL) new_indels += len;
      if (bam_cigar_type(op) & 1) q += len;
    }

    int64_t mismatches = md ? count_md_mismatches(md, gcigar, keep) : -1;
    if (mismatches < 0) {  // no usable MD
      if (!nm_tag) return;
      int64_t genome_indels = 0;
      for (uint32_t c : gcigar) {
        int op = bam_cigar_op(c);
        if (op == BAM_CINS || op == BAM_CDEL) genome_indels += bam_cigar_oplen(c);
      }
      mismatches = std::max<int64_t>(0, genome_nm - genome_indels);
    }

    int32_t nm = (int32_t)(mismatches + align.rescue_mm + new_indels);
    if (nm_tag) bam_aux_del(b, nm_tag);
    bam_aux_append(b, "NM", 'i', 4, (uint8_t*)&nm);
  }

  // Aux tags kept in the output. Everything else describes the genome
  // alignment (aligner scores, mate info, MD, ...) and is dropped. bramble's
  // own tags are set or replaced separately. Keep in sync with KEEP_TAGS in
  // bramble-cli/src/pipeline.rs.
  static const char* const KEEP_TAGS[] = {
    "NH", "HI", "NM", "AS",                      // bramble's own
    "RG",                                        // read group
    "BC", "QT", "CB", "CR", "CY", "MI", "OX",    // barcodes and UMIs
    "BZ", "RX", "QX", "UB", "UR", "UY", "BX",
    "MM", "ML", "MN", "Mm", "Ml",                // base modifications
  };

  // Remove every aux tag not in KEEP_TAGS
  void keep_only_known_tags(bam1_t* b) {
    uint8_t* a = bam_aux_first(b);
    while (a) {
      const char* t = bam_aux_tag(a);
      bool keep = false;
      for (const char* k : KEEP_TAGS) {
        if (t[0] == k[0] && t[1] == k[1]) { keep = true; break; }
      }
      a = keep ? bam_aux_next(b, a) : bam_aux_remove(b, a);
    }
  }

  // Put bramble's own tags first, in a fixed order (NH HI AS NM, as STAR
  // writes them), then the kept input tags in their input order. Same order as
  // the Rust version (pipeline.rs). Call after all tags are set.
  void order_tags(bam1_t* b) {
    static const char* const FIRST[] = {"NH", "HI", "AS", "NM"};
    uint8_t* aux = bam_get_aux(b);
    size_t len = bam_get_l_aux(b);
    if (len == 0) return;
    // [start, end) of each tag, tag name included
    thread_local std::vector<std::pair<size_t, size_t>> spans;
    spans.clear();
    uint8_t* a = bam_aux_first(b);
    while (a) {
      uint8_t* next = bam_aux_next(b, a);
      size_t start = (a - 2) - aux;
      size_t end = next ? (size_t)((next - 2) - aux) : len;
      spans.emplace_back(start, end);
      a = next;
    }
    thread_local std::vector<uint8_t> out;
    out.clear();
    auto is_first = [&](size_t start) {
      for (const char* t : FIRST)
        if (aux[start] == t[0] && aux[start + 1] == t[1]) return true;
      return false;
    };
    for (const char* t : FIRST)
      for (auto [s, e] : spans)
        if (aux[s] == t[0] && aux[s + 1] == t[1]) out.insert(out.end(), aux + s, aux + e);
    for (auto [s, e] : spans)
      if (!is_first(s)) out.insert(out.end(), aux + s, aux + e);
    if (out.size() == len) std::memcpy(aux, out.data(), len);
  }

  // Set mate information for a read
  void set_mate_info(bam1_t* b, BamInfo* this_pair, bool first_read) {

    // Handle unpaired reads
    if (!this_pair->is_paired) {
      // An orphan of a paired read stays paired (0x1, read1/read2 kept) with
      // its mate unmapped (0x8): the standard SAM orphan, as in the Rust
      // version. Single-end reads have no 0x1 and stay as they are.
      b->core.flag &= ~(BAM_FPROPER_PAIR | BAM_FMREVERSE);
      if (b->core.flag & BAM_FPAIRED) b->core.flag |= BAM_FMUNMAP;
      
      // Set mate reference to unmapped
      b->core.mtid = -1;  // -1 indicates unmapped mate
      b->core.mpos = -1;  // -1 indicates unmapped mate position
      b->core.isize = 0;  // No insert size for unpaired reads
      
      return;
    }
      
    // Set basic paired flags
    b->core.flag |= BAM_FPAIRED;
    
    // Mate-reverse flag = the mate record's output orientation: its genome
    // orientation, flipped when its transcript is on '-' (reverse_complement_bam).
    // Clear first: the genome alignment's 0x20 must not carry over.
    const auto &mate_out = first_read ? this_pair->read2 : this_pair->read1;
    char mate_strand = first_read ? this_pair->m_align.strand : this_pair->r_align.strand;
    bool mate_genome_rev = mate_out && mate_out->brec &&
      (mate_out->brec->get_b()->core.flag & BAM_FREVERSE);
    b->core.flag &= ~BAM_FMREVERSE;
    if (mate_genome_rev != (mate_strand == '-')) b->core.flag |= BAM_FMREVERSE;
    
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