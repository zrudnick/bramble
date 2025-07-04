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

uint32_t* get_new_cigar_array(uint32_t* cigar, uint32_t n_cigar, uint32_t* new_n_cigar) {
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
uint8_t* copy_cigar_memory(bam1_t* b, uint32_t new_n_cigar, uint32_t n_cigar, uint32_t new_m_data, 
    uint32_t* new_cigar, int l_qname, int l_qseq, int l_aux) {
    
    // Reallocate BAM data
    if (new_m_data > b->m_data) {
        b->m_data = new_m_data;
        b->data = (uint8_t*)realloc(b->data, b->m_data);
        if (!b->data) {
            GError("Failed to realloc bam data");
        }
    }
    uint8_t* data = b->data;

    // Copy new cigar replacing old cigar
    memcpy(data + l_qname, new_cigar, new_n_cigar * 4);

    int cigar_size = n_cigar * 4;
    int new_cigar_size = new_n_cigar * 4;
    int seq_offset = l_qname + cigar_size;
    int new_seq_offset = l_qname + new_cigar_size;
    int seq_size = (l_qseq + 1) / 2;
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
    int l_aux = b->l_data - (l_qname + n_cigar * 4 + (l_qseq + 1) / 2 + l_qseq);
    
    // Calculate new total data size
    uint32_t new_m_data = l_qname + new_n_cigar * 4 + (l_qseq + 1) / 2 + l_qseq + l_aux;

    copy_cigar_memory(b, new_n_cigar, n_cigar, new_m_data, new_cigar, l_qname, l_qseq, l_aux);
    b->core.n_cigar = new_n_cigar;

    free(new_cigar);
}

// Set mate information for a read
MateInfo* set_mate_info(BamIO* io, const std::string& transcript_name, tid_t transcript_id, bam1_t* b, int tid, 
    ReadInfo* read_info, g2tTree* g2t) {
    
    auto mate_info = read_info->mate_info[transcript_id];
    if (!mate_info || !mate_info->valid_pair) {
        // No mate information available or invalid pair
        b->core.mtid = -1;
        b->core.mpos = -1;
        b->core.isize = 0;
        // Don't set BAM_FPAIRED if no valid mate
        return nullptr;
    }
    
    // Set basic paired flags
    b->core.flag |= BAM_FPAIRED;
    
    if (mate_info->same_transcript) {
        // Both mates map to same transcript - use "=" for mate reference
        b->core.mtid = tid;  // Same as current read's tid
        b->core.mpos = mate_info->match_pos;
        b->core.flag |= BAM_FPROPER_PAIR;
        
        // Calculate insert size for same transcript
        int32_t isize = 0;
        if (b->core.pos < b->core.mpos) {
            // This read comes first
            isize = ((b->core.mpos + mate_info->mate_size) - b->core.pos);
        } else {
            // Mate comes first  
            isize = ((b->core.pos + read_info->read_size) - b->core.mpos);
            isize = -isize;  // Negative for second read in pair
        }
        b->core.isize = isize;
        
    } else {
        // Mates map to different transcripts

        // Get mate's transcript ID in BAM header
        auto& mate_transcript_name = g2t->get_tid_name(mate_info->transcript_id);
        int mate_tid = io->get_tid(mate_transcript_name.c_str());
        if (mate_tid < 0) {
            // Mate transcript not found in header - treat as unmapped mate
            // This shouldn't happen
            b->core.mtid = -1;
            b->core.mpos = -1;
            b->core.isize = 0;
            b->core.flag |= BAM_FMUNMAP;
            return nullptr;
        }
        
        // Set mate reference and position for different transcript
        b->core.mtid = mate_tid;
        b->core.mpos = mate_info->match_pos;
        b->core.isize = 0;
        b->core.flag |= BAM_FPROPER_PAIR;
    }
    
    // Set mate strand information
    if (read_info->mate_is_reverse) {
        b->core.flag |= BAM_FMREVERSE;
    }
     
    return mate_info;
}

// Add new NH tag (after removing current one)
void set_nh_tag(bam1_t* b, uint nh_i) {
    uint8_t* nh = bam_aux_get(b, "NH");
    if (nh) bam_aux_del(b, nh);
    bam_aux_append(b, "NH", 'i', sizeof(int32_t), (uint8_t*)&nh_i);
}

// Write bundle reads to BAM file
void write_to_bam(BamIO* io, std::map<tid_t, ReadInfo*>& bam_info, g2tTree* g2t) {
    GSamRecord new_brec;

    for (const auto& pair : bam_info) {
        auto read_info = std::get<1>(pair);
        if (!read_info || !read_info->valid_read) continue;

        // Get BAM and CIGAR string
        bam1_t* read_b = read_info->brec->get_b();
        uint32_t* cigar = bam_get_cigar(read_b);
        uint32_t n_cigar = read_b->core.n_cigar;

        // Remove any introns
        bool read_has_introns = has_introns(cigar, n_cigar);
        if (read_has_introns) update_cigar(read_b, cigar, n_cigar);

        // Set the NH (number of hits) tag
        auto nh_i = read_info->nh_i;
        set_nh_tag(read_b, nh_i);

        // Create new BAM record
        new_brec = *(read_info->brec);
        bam1_t* b = new_brec.get_b();

        for (const auto& match : read_info->matches) {
            tid_t transcript_id = std::get<0>(match);
            auto& transcript_name = g2t->get_tid_name(transcript_id);
            int tid = io->get_tid(transcript_name.c_str());
            uint match_pos = std::get<1>(match);
                
            // Set new coordinates
            b->core.tid = tid;
            b->core.pos = match_pos;

            // TODO: choose matches to be secondary

            // Set mate information if available
            auto mate_info = set_mate_info(io, transcript_name, transcript_id, b, tid, read_info, g2t);

            // Write the read match to the BAM
            io->write(&new_brec);

            // TODO: write mate at same time
            // Salmon needs to see paired reads next to each other
        }
    }  
}
