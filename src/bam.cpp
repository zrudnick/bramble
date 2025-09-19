
#include <set>
#include <unordered_set>
#include <string>
#include <vector>
#include <algorithm>
#include <memory>
#include <cstdlib>

#include "bramble.h"
#include "g2t.h"
#include "evaluate.h"
#include "bam.h"
#include "htslib/sam.h"

#ifndef NOTHREADS
#include "GThreads.h"
#endif

extern bool LONG_READS;
extern GFastMutex bam_io_mutex;  // Protects BAM io

namespace bramble {

// Add or merge a CIGAR operation
static uint32_t add_or_merge_op(uint32_t* new_cigar, uint32_t j, uint8_t op, uint32_t len) {
  // CIGAR array is not empty
  if (j > 0) {
    uint8_t prev_op = new_cigar[j-1] & BAM_CIGAR_MASK;
    if (prev_op == op) {
      uint32_t prev_len = new_cigar[j-1] >> BAM_CIGAR_SHIFT;
      new_cigar[j-1] = ((prev_len + len) << BAM_CIGAR_SHIFT) | op;
      return j;
    }
  }

  // CIGAR array is empty
  new_cigar[j] = (len << BAM_CIGAR_SHIFT) | op;
  return j + 1;
}

// Get new CIGAR array (for new soft clips)
uint32_t* get_new_cigar_soft_clips(uint32_t* cigar, uint32_t n_cigar, 
                              uint32_t* new_n_cigar, CigarMem& mem,
                              uint32_t soft_clip_front, uint32_t soft_clip_back) {

  uint32_t* new_cigar = mem.get_mem(n_cigar + 2); // allocate space for added soft clips
  uint32_t j = 0;
  
  uint32_t total_soft_clip_front = soft_clip_front;
  uint32_t total_soft_clip_back = soft_clip_back;
  
  // Check for existing front soft clip
  uint32_t start_i = 0;
  if ((cigar[0] & BAM_CIGAR_MASK) == BAM_CSOFT_CLIP) {
    total_soft_clip_front += (cigar[0] >> BAM_CIGAR_SHIFT);
    start_i = 1;
  }
  
  // Check for existing back soft clip
  uint32_t end_i = n_cigar;
  if ((cigar[n_cigar-1] & BAM_CIGAR_MASK) == BAM_CSOFT_CLIP) {
    total_soft_clip_back += (cigar[n_cigar-1] >> BAM_CIGAR_SHIFT);
    end_i = n_cigar - 1;
  }

  // Add front soft clip
  if (total_soft_clip_front > 0) {
    j = add_or_merge_op(new_cigar, j, BAM_CSOFT_CLIP, total_soft_clip_front);
  }
  
  // Find where front clipping ends
  while (start_i < end_i && soft_clip_front > 0) {
    uint8_t op = cigar[start_i] & BAM_CIGAR_MASK;
    if (op == BAM_CREF_SKIP) {
      continue;
    }
    
    uint32_t len = cigar[start_i] >> BAM_CIGAR_SHIFT;
    if (soft_clip_front >= len) {
      soft_clip_front -= len;
      start_i++;
    } else {
      break; // op will be partially preserved
    }
  }
  
  // Find where the back clipping starts
  int32_t test_i = end_i - 1;
  while (test_i >= (int32_t)start_i && soft_clip_back > 0) {
    uint8_t op = cigar[test_i] & BAM_CIGAR_MASK;
    if (op == BAM_CREF_SKIP) {
      test_i--;
      end_i--;
      continue;
    }
    
    uint32_t len = cigar[test_i] >> BAM_CIGAR_SHIFT;
    if (soft_clip_back >= len) {
      soft_clip_back -= len;
      end_i--;
      test_i--;
    } else {
      break; // op will be partially preserved
    }
  }
  
  // Process the remaining operations
  for (uint32_t i = start_i; i < end_i; i++) {
    uint8_t op = cigar[i] & BAM_CIGAR_MASK;
    
    // Skip N operations
    if (op == BAM_CREF_SKIP) {
      continue; 
    }
    
    uint32_t len = cigar[i] >> BAM_CIGAR_SHIFT;
    
    // Partial front clipping
    if (i == start_i && soft_clip_front > 0) {
      len -= soft_clip_front;
    }
    
    // Back clipping on this operation
    if (i == (end_i - 1) && soft_clip_back > 0) {
      len -= soft_clip_back;
    }
    
    if (len > 0) {
      j = add_or_merge_op(new_cigar, j, op, len);
    }
  }

  // Add back soft clip
  if (total_soft_clip_back > 0) {
    j = add_or_merge_op(new_cigar, j, BAM_CSOFT_CLIP, total_soft_clip_back);
  }

  *new_n_cigar = j;
  return new_cigar;
}

// Get new CIGAR array (no new soft clips)
uint32_t* get_new_cigar(uint32_t* cigar, uint32_t n_cigar, 
                        uint32_t* new_n_cigar, CigarMem& mem) {
  
  uint32_t* new_cigar = mem.get_mem(n_cigar); 
  uint32_t j = 0;
  
  for (uint32_t i = 0; i < n_cigar; ++i) {
    uint8_t op = cigar[i] & BAM_CIGAR_MASK;
    uint32_t len = cigar[i] >> BAM_CIGAR_SHIFT;
    
    if (op != BAM_CREF_SKIP) {  // Skip N operations
      // Check if we can merge with the previous operation
      if (j > 0) {
        uint8_t prev_op = new_cigar[j-1] & BAM_CIGAR_MASK;
        uint32_t prev_len = new_cigar[j-1] >> BAM_CIGAR_SHIFT;
        
        // Merge with previous operation of same type
        if (prev_op == op) {
          new_cigar[j-1] = ((prev_len + len) << BAM_CIGAR_SHIFT) | op;

        // Different operation type, add new entry
        } else {
          new_cigar[j] = (len << BAM_CIGAR_SHIFT) | op;
          j++;
        }
      
      // First operation, just add it
      } else {
        new_cigar[j] = (len << BAM_CIGAR_SHIFT) | op;
        j++;
      }
    }
    // N operations are skipped (not copied)
  }

  *new_n_cigar = j;
  return new_cigar;
}

// Copy CIGAR memory for intron removal
uint8_t* copy_cigar_memory(bam1_t* b, uint32_t new_n_cigar, uint32_t old_n_cigar,
                           uint32_t new_m_data, const uint32_t* new_cigar,
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
                  CigarMem& mem, uint32_t soft_clip_front,
                  uint32_t soft_clip_back) {

  uint32_t* new_cigar = nullptr;
  uint32_t new_n_cigar = 0;

  // Build a new cigar (with or without soft clips)
  if (soft_clip_front || soft_clip_back) {
      new_cigar = get_new_cigar_soft_clips(cigar, n_cigar, &new_n_cigar,
                                           mem, soft_clip_front, soft_clip_back);
  } else {
      new_cigar = get_new_cigar(cigar, n_cigar, &new_n_cigar, mem);
  }

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

void write_to_bam(BamIO *io, std::unordered_map<bam_id_t, BamInfo *>& bam_info) {
  std::unordered_set<read_id_t> seen;
  CigarMem cigar_mem; 

  // Print information for both mates (paired reads) or single read
  for (const auto& [id, this_pair] : bam_info) {
    if (!this_pair || !this_pair->valid_pair) continue;

    auto write_read = [&](ReadOut* read, bool is_first) {
      if (!read || !read->brec) return;
      
      bam1_t* b = read->brec->get_b();

      if (seen.insert(read->index).second) {
        
        uint32_t* cigar = bam_get_cigar(b);
        uint32_t n_cigar = b->core.n_cigar;

        update_cigar(b, cigar, n_cigar, cigar_mem, 
                     read->soft_clip_front, read->soft_clip_back);

        // Set NH and XS tags
        set_nh_tag(b, read->nh);
        set_xs_tag(b, read->strand);

        // Set FLAG to indicate whether reverse
        if (read->is_reverse) b->core.flag |= BAM_FREVERSE;

      }

      // Set new coordinates
      if (is_first) {
        b->core.tid = (int32_t)this_pair->r_tid;
        b->core.pos = (int32_t)this_pair->r_pos;
      } else {
        b->core.tid = (int32_t)this_pair->m_tid;
        b->core.pos = (int32_t)this_pair->m_pos;
      }

      // Set mate information if available
      if (this_pair->is_paired)
          set_mate_info(b, this_pair, is_first);

#ifndef NOTHREADS
      bam_io_mutex.lock();
#endif

      io->write(read->brec);

#ifndef NOTHREADS
      bam_io_mutex.unlock();
#endif
    };

    // Always process first read
    write_read(this_pair->read1, true);     // first_read = true

    // Process second read if paired
    if (this_pair->is_paired) {
      write_read(this_pair->read2, false);  // first_read = false
    }

  }
}

} // namespace bramble