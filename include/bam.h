
#pragma once

using namespace bramble;

namespace bramble {

  struct BamInfo;
  struct BamIO;
  struct Cigar;

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

  // -------- function definitions

  uint32_t* get_new_cigar(uint32_t* real_cigar, uint32_t n_real_cigar,
                        const Cigar& ideal_cigar, uint32_t* new_n_cigar, 
                        CigarMem& mem);

  uint8_t* copy_cigar_memory(bam1_t* b, uint32_t new_n_cigar,
                            uint32_t old_n_cigar, uint32_t new_m_data,
                            const uint32_t* new_cigar, int l_qname,
                            int l_qseq, int l_aux);

  bool update_cigar(bam1_t* b, uint32_t* cigar, uint32_t n_cigar,
                    CigarMem& mem, const Cigar& ideal_cigar);

  void set_mate_info(bam1_t* b, BamInfo* this_pair, bool first_read);

  void set_nh_tag(bam1_t* b, int32_t nh_i);

  void set_xs_tag(bam1_t* b, char xs_a);

  void remove_extra_tags(bam1_t* b);

}