// process_reads.h

#ifndef BAM_H
#define BAM_H

#include <memory>

#include "bramble.h"

using tid_t = uint32_t;
using read_id_t = uint32_t;

bool has_introns(uint32_t* cigar, uint32_t n_cigar);

uint32_t* get_new_cigar_array(uint32_t* cigar, uint32_t n_cigar, uint32_t* new_n_cigar);

uint8_t* copy_cigar_memory(bam1_t* b, uint32_t new_n_cigar, uint32_t n_cigar, uint32_t new_m_data, uint32_t* new_cigar,
                           int l_qname, int l_qseq, int l_aux);

void update_cigar(bam1_t* b, uint32_t* cigar, uint32_t n_cigar);

MateInfo* set_mate_info(BamIO* io, tid_t tid, bam1_t* b, ReadInfo* read_info, g2tTree* gt2);

void set_nh_tag(bam1_t* b, uint nh_i);

void write_to_bam(BamIO* io, std::map<tid_t, ReadInfo*>& bam_info, g2tTree* gt2);

#endif
