// process_reads.h

#ifndef BAM_H
#define BAM_H

#include <memory>

#include "bramble.h"

using tid_t = uint32_t;
using read_id_t = uint32_t;

uint32_t* get_new_cigar_array(uint32_t* cigar, uint32_t n_cigar, uint32_t* new_n_cigar);

uint8_t* copy_cigar_memory(bam1_t* b, uint32_t new_n_cigar, uint32_t n_cigar, uint32_t new_m_data, 
                           uint32_t* new_cigar, int l_qname, int l_qseq, int l_aux);

void update_cigar(bam1_t* b, uint32_t* cigar, uint32_t n_cigar);

void set_mate_info(bam1_t* b, BamInfo* this_pair, bool first_read);

void set_nh_tag(bam1_t* b, uint nh_i);

void write_to_bam(BamIO* io, std::map<bam_id_t, BamInfo*>& bam_info, g2tTree* g2t);

#endif
