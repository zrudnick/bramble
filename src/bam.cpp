// bam.cpp

#include <set>
#include <string>
#include <vector>
#include <algorithm>
#include <memory>

#include "bramble.h"
#include "htslib/sam.h"

// Check if CIGAR string contains introns
bool has_introns(uint32_t* cigar, uint32_t n_cigar) {
    for (uint32_t i = 0; i < n_cigar; i++) {
        if ((cigar[i] & BAM_CIGAR_MASK) == BAM_CREF_SKIP) {
            return true;
        }
    }
    return false;
}

uint32_t* get_new_cigar_array(uint32_t* cigar, uint32_t n_cigar, 
                              uint32_t* new_n_cigar) {
  uint32_t* new_cigar = (uint32_t*)malloc(n_cigar * sizeof(uint32_t));    // max new_n_cigar = n_cigar
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
    // N operations are simply skipped (not copied)
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
void update_cigar(bam1_t* b, uint32_t* cigar, uint32_t n_cigar) {
  // Create new CIGAR array without N operations
  uint32_t new_n_cigar;
  uint32_t* new_cigar = get_new_cigar_array(cigar, n_cigar, &new_n_cigar);

  // Calculate new data layout sizes
  int l_qname = b->core.l_qname;
  int l_qseq = b->core.l_qseq;
  int l_aux = b->l_data - (l_qname + (n_cigar << 2) + ((l_qseq + 1) >> 1) + l_qseq);
  
  // Calculate new total data size
  uint32_t new_m_data = l_qname + (new_n_cigar << 2) + ((l_qseq + 1) >> 1) + l_qseq + l_aux;

  copy_cigar_memory(b, new_n_cigar, n_cigar, new_m_data, new_cigar, l_qname, l_qseq, l_aux);
  b->core.n_cigar = new_n_cigar;

  free(new_cigar);
}

// Set mate information for a read
void set_mate_info(bam1_t* b, BamInfo* this_pair, bool first_read) {
    
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
void write_to_bam(BamIO* io, std::map<bam_id_t, BamInfo*>& bam_info, g2tTree* g2t) {
  GSamRecord new_brec;
  std::set<read_id_t> seen;

  // Each pair now contains info for both mates
  // So they are printed at the same time
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

      // Remove any introns
      if (has_introns(cigar, n_cigar)) {
        update_cigar(read_b, cigar, n_cigar); // TODO: is there a way to speed this up?
      }

      // Set the NH (number of hits) tag
      nh_i = this_pair->nh;
      set_nh_tag(read_b, nh_i);

      seen.insert(this_pair->read_index);
    }

    // Create new BAM record
    new_brec = *(this_pair->brec);
    bam1_t* b = new_brec.get_b();

    // Set new coordinates
    b->core.tid = this_pair->tid;
    b->core.pos = this_pair->pos;

    // Set mate information if available
    set_mate_info(b, this_pair, true);

    io->write(&new_brec);

    // 2. Process second mate
    if (seen.find(this_pair->mate_index) == seen.end()) {
      // Get BAM and CIGAR string
      read_b = this_pair->mate_brec->get_b();

      cigar = bam_get_cigar(read_b);
      n_cigar = read_b->core.n_cigar;

      // Remove any introns
      if (has_introns(cigar, n_cigar)) {
        update_cigar(read_b, cigar, n_cigar);
      }
  
      // Set the NH (number of hits) tag
      nh_i = this_pair->mate_nh;
      set_nh_tag(read_b, nh_i);

      seen.insert(this_pair->mate_index);
    }

    // Create new BAM record
    new_brec = *(this_pair->mate_brec);
    b = new_brec.get_b();

    // Set new coordinates
    b->core.tid = this_pair->mate_tid;
    b->core.pos = this_pair->mate_pos;

    // Set mate information if available
    set_mate_info(b, this_pair, false);

    io->write(&new_brec);
  }

  seen.clear();
}
