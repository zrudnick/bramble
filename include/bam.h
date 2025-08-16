// process_reads.h

#ifndef BAM_H
#define BAM_H

#include <memory>

#include "bramble.h"

using tid_t = uint32_t;
using read_id_t = uint32_t;

struct CigarMem {
  uint32_t* ptr = nullptr;
  size_t len = 0;

  ~CigarMem() {
    if (len > 0) {
      if (ptr != nullptr) {
        free(ptr);
        ptr = nullptr;
      }
    }
  }
  void clear() { 
    if (len > 0) {
      memset(reinterpret_cast<char*>(ptr), 0, sizeof(uint32_t) * len);
    }
  }

  uint32_t* get_mem(uint32_t size) {
    if (len == 0) {
      uint32_t claim = size > 0 ? size : 1;
      ptr = reinterpret_cast<uint32_t*>(malloc(sizeof(uint32_t) * size));
      if (ptr == nullptr) { std::exit(EXIT_FAILURE); }
      len = claim;
    } else if (len < size) { 
      ptr = reinterpret_cast<uint32_t*>(realloc(ptr, sizeof(uint32_t) * size));
      if (ptr == nullptr) { std::exit(EXIT_FAILURE); }
      len = size;
    } 
    return ptr;
  }
};

uint32_t* get_new_cigar_array(uint32_t* cigar, uint32_t n_cigar, uint32_t* new_n_cigar);

uint8_t* copy_cigar_memory(bam1_t* b, uint32_t new_n_cigar, uint32_t n_cigar, uint32_t new_m_data, 
                           uint32_t* new_cigar, int l_qname, int l_qseq, int l_aux);

bool update_cigar(bam1_t* b, uint32_t* cigar, uint32_t n_cigar, CigarMem& mem, 
                  uint32_t soft_clip_front, uint32_t soft_clip_back);

void set_mate_info(bam1_t* b, BamInfo* this_pair, bool first_read);

void set_nh_tag(bam1_t* b, uint nh_i);

void write_to_bam(BamIO* io, std::map<bam_id_t, BamInfo*>& bam_info, g2tTree* g2t);

#endif
