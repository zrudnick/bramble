// bam.cpp

#include <set>
#include <unordered_set>
#include <string>
#include <vector>
#include <algorithm>
#include <memory>
#include <cstdlib>

#include "bramble.h"
#include "bam.h"
#include "htslib/sam.h"

#ifndef NOTHREADS
#include "GThreads.h"
#endif

GFastMutex bam_io_mutex;  // Protects writes to BAM output

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

// Get new CIGAR array (if there is a new soft clip to add)
uint32_t* get_new_cigar_array_soft_clips(uint32_t* cigar, uint32_t n_cigar, 
                              uint32_t* new_n_cigar, CigarMem& mem,
                              uint32_t soft_clip_front, uint32_t soft_clip_back,
                              bool &discard) {
    
  uint32_t* new_cigar = mem.get_mem(n_cigar);
  uint32_t j = 0;

  uint32_t total_soft_clip_front = soft_clip_front;
  uint32_t total_soft_clip_back = soft_clip_back;
  
  // Check for existing front soft clip
  uint32_t start_i = 0;
  if ((cigar[0] & BAM_CIGAR_MASK) == BAM_CSOFT_CLIP) {
      total_soft_clip_front += (cigar[0] >> BAM_CIGAR_SHIFT); // # existing soft clip at front
      start_i = 1;
  }
  
  // Check for existing back soft clip
  uint32_t end_i = n_cigar;
  if ((cigar[n_cigar-1] & BAM_CIGAR_MASK) == BAM_CSOFT_CLIP) {
    total_soft_clip_back += (cigar[n_cigar-1] >> BAM_CIGAR_SHIFT); // # existing soft clip at back
    end_i = n_cigar - 1;
  }

  // Check if soft clip will be adjacent to insertion/deletion at the beginning
  if (soft_clip_front > 0 && start_i < end_i) {
      uint8_t first_op = cigar[start_i] & BAM_CIGAR_MASK;
      if (first_op == BAM_CDEL || first_op == BAM_CINS) {
          discard = true;
          *new_n_cigar = 0;
          return new_cigar;
      }
  }
  
  // Check if soft clip will be adjacent to insertion/deletion at the end
  if (soft_clip_back > 0 && start_i < end_i) {
      uint8_t last_op = cigar[end_i-1] & BAM_CIGAR_MASK;
      if (last_op == BAM_CDEL || last_op == BAM_CINS) {
          discard = true;
          *new_n_cigar = 0;
          return new_cigar;
      }
  }

  // Add front soft clip
  if (total_soft_clip_front > 0) {
    j = add_or_merge_op(new_cigar, j, BAM_CSOFT_CLIP, total_soft_clip_front);
  }
  
  // Process middle operations, skipping N and adjusting inside ops
  uint32_t i = start_i;
  uint32_t to_subtract_front = soft_clip_front;
  uint32_t to_subtract_back = soft_clip_back;

  while (i < end_i) {
    
    // Get current op in CIGAR
    uint8_t op = cigar[i] & BAM_CIGAR_MASK;

    if (op == BAM_CREF_SKIP) {
      i++;
      continue;   // don't add N operations
    }

    // Get current length of op
    uint32_t len = cigar[i] >> BAM_CIGAR_SHIFT;
    
    // Adjust first operation (subtract front soft clip)
    if (to_subtract_front > 0) {
      // Discard read if deletion/insertion intersects soft clip
      if (op == BAM_CDEL || op == BAM_CINS) {
        discard = true;
        *new_n_cigar = 0;
        return new_cigar;
      }

      if (to_subtract_front >= len) {
        // Soft clip covers this op completely, don't add
        to_subtract_front -= len;
        i++;
      } else {
        // Soft clip reduces this op's size
        len -= to_subtract_front;
        to_subtract_front = 0;  // Clear the remaining front clip
        j = add_or_merge_op(new_cigar, j, op, len);
        i++;
      }  
    } else if (i < end_i - 1 || to_subtract_back == 0) {
      // Normal case: add the operation as-is
      j = add_or_merge_op(new_cigar, j, op, len);
      i++;
    } else {
      // Last operation and we need to subtract back soft clip
      // Discard read if deletion/insertion intersects soft clip
      if (op == BAM_CDEL || op == BAM_CINS) {
        discard = true;
        *new_n_cigar = 0;
        return new_cigar;
      }
      
      if (to_subtract_back >= len) {
        // Soft clip covers this op completely, don't add
        // Move to previous operation
        to_subtract_back -= len;
        end_i--;
        if (i >= end_i) {
          // We've gone past all operations, handle remaining soft clip
          break;
        }
        // Don't increment i, reprocess current position
      } else {
        // Soft clip reduces this op's size
        len -= to_subtract_back;
        to_subtract_back = 0;  // Clear the remaining back clip
        j = add_or_merge_op(new_cigar, j, op, len);
        i++;
      } 
    }
  }

  // Handle any remaining back soft clip that couldn't be subtracted
  // This happens when the soft clip is larger than available operations
  if (to_subtract_back > 0) {
      // Need to go back through already added operations
      uint32_t remaining_clip = to_subtract_back;
      int32_t back_idx = j - 1;
      
      while (back_idx >= 0 && remaining_clip > 0) {
          uint8_t back_op = new_cigar[back_idx] & BAM_CIGAR_MASK;
          uint32_t back_len = new_cigar[back_idx] >> BAM_CIGAR_SHIFT;
          
          // Check if this operation is an insertion or deletion
          if (back_op == BAM_CDEL || back_op == BAM_CINS) {
              discard = true;
              *new_n_cigar = 0;
              return new_cigar;
          }
          
          if (remaining_clip >= back_len) {
              // Remove this operation entirely
              remaining_clip -= back_len;
              j--;
              back_idx--;
          } else {
              // Reduce this operation's length
              uint32_t new_len = back_len - remaining_clip;
              new_cigar[back_idx] = (new_len << BAM_CIGAR_SHIFT) | back_op;
              remaining_clip = 0;
          }
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
uint32_t* get_new_cigar_array(uint32_t* cigar, uint32_t n_cigar, 
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
uint8_t* copy_cigar_memory(bam1_t* b, uint32_t new_n_cigar, 
                           uint32_t n_cigar, uint32_t new_m_data, 
                           uint32_t* new_cigar, int l_qname, 
                           int l_qseq, int l_aux) {
    
  // Reallocate BAM data
  if (new_m_data > b->m_data) {
      b->m_data = new_m_data;
      b->data = (uint8_t*)realloc(b->data, b->m_data);
  }
  uint8_t* data = b->data;

  // Copy new cigar replacing old cigar
  memcpy(data + l_qname, new_cigar, new_n_cigar << 2);

  int seq_offset = l_qname + (n_cigar << 2);
  int new_seq_offset = l_qname + (new_n_cigar << 2);
  int seq_size = (l_qseq + 1) >> 1;
  int qual_size = l_qseq;

  memmove(data + new_seq_offset, data + seq_offset, seq_size + qual_size + l_aux);
  b->l_data = new_m_data;

  return data;
}

// Update CIGAR to remove introns
bool update_cigar(bam1_t* b, uint32_t* cigar, uint32_t n_cigar, CigarMem& mem, 
                  uint32_t soft_clip_front, uint32_t soft_clip_back) {
  
  // Create new CIGAR array without N operations
  uint32_t new_n_cigar;

  uint32_t* new_cigar;
  if (soft_clip_front || soft_clip_back) {
    bool discard = false;
    new_cigar = get_new_cigar_array_soft_clips(cigar, n_cigar, &new_n_cigar, mem, 
                          soft_clip_front, soft_clip_back, discard);
    if (discard) return false;
  } else {
    new_cigar = get_new_cigar_array(cigar, n_cigar, &new_n_cigar, mem);
  }

  // Calculate new data layout sizes
  int l_qname = b->core.l_qname;
  int l_qseq = b->core.l_qseq;
  int l_aux = b->l_data - (l_qname + (n_cigar << 2) + ((l_qseq + 1) >> 1) + l_qseq);
  
  // Calculate new total data size
  uint32_t new_m_data = l_qname + (new_n_cigar << 2) + ((l_qseq + 1) >> 1) + l_qseq + l_aux;

  copy_cigar_memory(b, new_n_cigar, n_cigar, new_m_data, new_cigar, l_qname, l_qseq, l_aux);
  b->core.n_cigar = new_n_cigar;

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
    tid = this_pair->mate_tid;
    pos = this_pair->mate_pos;
  } else {
    tid = this_pair->tid;
    pos = this_pair->pos;
  }
  
  // Set basic paired flags
  b->core.flag |= BAM_FPAIRED;
  
  // Mates map to same transcript
  if (this_pair->same_transcript) {
    // Both mates map to same transcript - use "=" for mate reference 
    b->core.mtid = tid; // Same as current read's tid
    b->core.mpos = pos;
    b->core.flag |= BAM_FPROPER_PAIR;
    
    // Calculate insert size for same transcript
    int32_t isize = 0;
    if (first_read) {
      // This read comes first
      isize = ((b->core.mpos + this_pair->mate_size) - b->core.pos);
    } else {
      // Mate comes first
      isize = -((b->core.pos + this_pair->mate_size) - b->core.mpos); // negative for second read in pair
    }
    b->core.isize = isize;
      
  // Mates map to different transcripts
  } else {
    // Set mate reference and position for different transcript
    b->core.mtid = tid;
    b->core.mpos = pos;
    b->core.isize = 0;
    b->core.flag |= BAM_FPROPER_PAIR;
  }
  
  // Set mate strand information
  if (first_read && this_pair->mate_is_reverse) {
    b->core.flag |= BAM_FMREVERSE;
  } else if (!first_read && this_pair->is_reverse) {
    b->core.flag |= BAM_FMREVERSE;
  }
}

// Add new NH tag (after removing current one)
void set_nh_tag(bam1_t* b, uint nh_i) {
  uint8_t* nh = bam_aux_get(b, "NH");
  if (nh) bam_aux_del(b, nh);
  bam_aux_append(b, "NH", 'i', sizeof(int32_t), (uint8_t*)&nh_i);
}

// Write bundle reads to BAM file
void write_to_bam(BamIO* io, std::map<bam_id_t, BamInfo*>& bam_info) {
  GSamRecord new_brec;
  std::unordered_set<read_id_t> seen;

  CigarMem cigar_mem; 

  // Print information for both mates (paired reads) or single read
  for (const auto& pair : bam_info) {
    auto this_pair = std::get<1>(pair);
    if (!this_pair || !this_pair->valid_pair) continue;

    bam1_t* read_b;
    uint32_t* cigar;
    uint32_t n_cigar;
    uint32_t nh_i;

    // 1. Process first mate
    if (seen.find(this_pair->read_index) == seen.end()) {
      // Get BAM and CIGAR string
      read_b = this_pair->brec->get_b();
      cigar = bam_get_cigar(read_b);
      n_cigar = read_b->core.n_cigar;
      
      uint32_t soft_clip_front = this_pair->soft_clip_front;
      uint32_t soft_clip_back = this_pair->soft_clip_back;

      if (update_cigar(read_b, cigar, n_cigar, cigar_mem, soft_clip_front, soft_clip_back)) {
        // Set the NH (number of hits) tag
        nh_i = this_pair->nh;
        set_nh_tag(read_b, nh_i);
      } else {
        this_pair->discard_read = true;
      }

      seen.insert(this_pair->read_index);
    }

    // Write valid read
    if (!this_pair->discard_read) {
      // Create new BAM record
      new_brec = *(this_pair->brec);
      bam1_t* b = new_brec.get_b();

      // Set new coordinates
      b->core.tid = this_pair->tid;
      b->core.pos = this_pair->pos;

      // Set mate information if available
      set_mate_info(b, this_pair, true);

      bam_io_mutex.lock();
      io->write(&new_brec);
      bam_io_mutex.unlock();
    }
    
    if (!this_pair->is_paired) continue;  // skip mate processing if unpaired

    // 2. Process second mate
    if (seen.find(this_pair->mate_index) == seen.end()) {
      // Get BAM and CIGAR string
      read_b = this_pair->mate_brec->get_b();

      cigar = bam_get_cigar(read_b);
      n_cigar = read_b->core.n_cigar;

      uint32_t soft_clip_front = this_pair->mate_soft_clip_front;
      uint32_t soft_clip_back = this_pair->mate_soft_clip_back;

      if (update_cigar(read_b, cigar, n_cigar, cigar_mem, soft_clip_front, soft_clip_back)) {

        // Set the NH (number of hits) tag
        nh_i = this_pair->mate_nh;
        set_nh_tag(read_b, nh_i);
      } else {
        this_pair->mate_discard_read = true;
      }

      seen.insert(this_pair->mate_index);
    }

    if (!this_pair->mate_discard_read) {
      // Create new BAM record
      new_brec = *(this_pair->mate_brec);
      bam1_t* b = new_brec.get_b();

      // Set new coordinates
      b->core.tid = this_pair->mate_tid;
      b->core.pos = this_pair->mate_pos;

      // Set mate information if available
      set_mate_info(b, this_pair, false);

      bam_io_mutex.lock();
      io->write(&new_brec);
      bam_io_mutex.unlock();
    }  

  }

  seen.clear();
}
