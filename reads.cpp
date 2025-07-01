// reads.cpp

#include "bramble.h"
#include "htslib/sam.h"

extern bool longreads;
extern uint junctionsupport; 	// anchor length for junction to be considered well supported

static GStr read_id("", 256);   // to prevent repeated reallocation for each parsed read
//not thread safe -- to only be used in processRead() as long as that's the unique producer

bool has_deletion_at_junction(GSamRecord* brec, int junction_index) {
    return (brec->juncsdel[junction_index].start || 
            brec->juncsdel[junction_index].end);
}

bool has_guide_match(GReadAlnData& alndata, int junction_index) {
    return (alndata.juncs.Count() && alndata.juncs[junction_index]->guide_match);
}

CJunction* create_junction(GSamRecord* brec, int exon_index, 
        char strand, GReadAlnData& alndata, GList<CJunction>& junction) {
    
    int junction_strand = strand;
    uint junction_start = brec->exons[exon_index - 1].end;
    uint junction_end = brec->exons[exon_index].start;
    
    // Handle long read deletions at junction boundaries
    if (longreads && has_deletion_at_junction(brec, exon_index - 1)) {
        if (!has_guide_match(alndata, exon_index - 1)) {
            junction_strand = 0;
            junction_start -= brec->juncsdel[exon_index - 1].start;
            junction_end += brec->juncsdel[exon_index - 1].end;
        }
    }
    
    CJunction* new_junction = junction.AddIfNew(
        new CJunction(junction_start, junction_end, junction_strand), true);
    
    if (alndata.juncs.Count() && new_junction) {
        new_junction->guide_match = alndata.juncs[exon_index - 1]->guide_match;
    }
    
    return new_junction;
}

void process_exons_and_junctions(GSamRecord* brec, CReadAln* readaln,
        GList<CJunction>& junction, GReadAlnData& alndata) {
    auto exons = brec->exons;
    for (int i = 0; i < exons.Count(); i++) {
        readaln->len += exons[i].len();
        
        if (i > 0) {
            // Ensure null junction exists
            if (junction.Count() == 0) {
        		CJunction* null_junction = new CJunction(0, 0, 0);
        		junction.Add(null_junction);
    		}
            
            // Create junction
            CJunction* new_junction = create_junction(brec, i, readaln->strand, alndata, junction);
            
            if (new_junction) {
                readaln->juncs.Add(new_junction);
            }
        }
        
        readaln->segs.Add(exons[i]);
    }
}

void update_bundle_end(uint bundle_end, BundleData& bundle, int read_end) {
    if (read_end > bundle_end) {
        bundle_end = read_end;
        bundle.end = bundle_end;
    }
}

float calculate_read_count(GSamRecord* brec, float unitig_cov) {
    
    float read_count = static_cast<float>(brec->tag_int("YC"));
    if (read_count == 0) read_count = 1;
    if (unitig_cov) read_count = unitig_cov;
    return read_count;
}

double adjust_mismatch_count(GSamRecord* brec, double nm) {
    if (!nm) {
        nm = static_cast<double>(brec->tag_int("nM"));
        if (brec->isPaired()) {
            nm /= 2;
        }
    }
    
    if (brec->clipL) nm++;
    if (brec->clipR) nm++;
    
    return nm;
}

bool process_matching_bases(char*& parser, uint& parsed_length) {
    unsigned int num_matches = 0;
    parseUInt(parser, num_matches);
    parsed_length += num_matches;
    return true;
}

void update_segment_tracking(CReadAln* read_aln, uint& segment_index, uint& read_length, int parsed_length) {
    while (segment_index < read_aln->segs.Count() && 
           read_length + static_cast<int>(read_aln->segs[segment_index].len()) < parsed_length) {
        read_length += read_aln->segs[segment_index].len();
        segment_index++;
    }
}

bool is_deletion_near_splice_site(CReadAln* read_aln, int segment_index, int parsed_length, int deletion_length) {
    // Check left splice site proximity
    if (segment_index > 0 && 
        parsed_length - read_aln->segs[segment_index].start < junctionsupport) {
        return true;
    }
    
    // Check right splice site proximity
    if (segment_index < read_aln->segs.Count() - 1 && 
        read_aln->segs[segment_index].end + 1 - parsed_length - deletion_length < junctionsupport) {
        return true;
    }
    
    return false;
}

bool process_deletion(CReadAln* read_aln, char*& c, uint& segment_index, 
                    uint& parsed_length, uint& read_length) {
    
    uint deletion_length = 0;
    char del_c = *(++c);
    
    // Count deletion length
    while (del_c >= 'A' && del_c <= 'Z') {
        deletion_length++;
        del_c = *(++c);
    }
    
    // Update segment tracking
    update_segment_tracking(read_aln, segment_index, read_length, parsed_length);
    
    if (segment_index == read_aln->segs.Count()) return false;
    
    // Check if deletion is too close to splice sites
    if (is_deletion_near_splice_site(read_aln, segment_index, parsed_length, deletion_length)) {
        return true;
    }
    
    parsed_length += deletion_length;
    return false;
}

bool is_mismatch_near_splice_site(CReadAln* read_aln, int segment_index, int parsed_length) {
    // Check left splice site proximity
    if (segment_index > 0 && 
        parsed_length - read_aln->segs[segment_index].start < junctionsupport) {
        return true;
    }
    
    // Check right splice site proximity
    if (segment_index < read_aln->segs.Count() - 1 && 
        read_aln->segs[segment_index].end - parsed_length < junctionsupport) {
        return true;
    }
    
    return false;
}

bool process_mismatch(CReadAln* read_aln, uint& segment_index, int parsed_length, uint& read_length) {
    
    update_segment_tracking(read_aln, segment_index, read_length, parsed_length);
    if (segment_index == read_aln->segs.Count()) return false;
    
    return is_mismatch_near_splice_site(read_aln, segment_index, parsed_length);
}

bool process_MD_string(CReadAln* read_aln, char* md_string, int ref_start) {
    char* c = md_string;
    uint segment_index = 0;
    uint parsed_length = ref_start;
    uint read_length = 0;
    
    while (*c != '\0') {
        if (*c >= '0' && *c <= '9') {
            if (process_matching_bases(c, parsed_length)) continue;
        }
        
        if (*c == '^') {
            if (process_deletion(read_aln, c, segment_index, parsed_length, read_length)) {
                return true;
            }
            continue;
        }
        
        if (*c >= 'A' && *c <= 'Z') {
            if (process_mismatch(read_aln, segment_index, parsed_length, read_length)) {
                return true;
            }
            parsed_length++;
        }
        
        c++;
    }
    
    return false;
}

bool is_insertion_near_splice_site(CReadAln* read_aln, int segment_index, int parsed_length) {
    // Check left splice site proximity
    if (segment_index > 0 && 
        parsed_length - read_aln->segs[segment_index].start < junctionsupport) {
        return true;
    }
    
    // Check right splice site proximity
    if (segment_index < read_aln->segs.Count() - 1 && 
        read_aln->segs[segment_index].end - parsed_length < junctionsupport) {
        return true;
    }
    
    return false;
}

bool process_insertion(CReadAln* read_aln, uint& segment_index, int parsed_length, uint& read_length) {
    update_segment_tracking(read_aln, segment_index, read_length, parsed_length);
    
    if (segment_index == read_aln->segs.Count()) {
        return false;
    }
    
    return is_insertion_near_splice_site(read_aln, segment_index, parsed_length);
}

bool process_cigar_string(CReadAln* read_aln, bam1_t* bam_record) {
    uint32_t* cigar = bam_get_cigar(bam_record);
    uint read_length = 0;
    uint parsed_length = 0;
    uint segment_index = 0;
    
    for (uint j = 0; j < bam_record->core.n_cigar; ++j) {
        int operation = bam_cigar_op(cigar[j]);
        
        if ((operation == BAM_CMATCH || operation == BAM_CEQUAL ||
             operation == BAM_CDIFF || operation == BAM_CDEL)) {
            parsed_length += bam_cigar_oplen(cigar[j]);
        } else if (operation == BAM_CINS) {
            if (process_insertion(read_aln, segment_index, parsed_length, read_length)) {
                return true;
            }
        }
    }
    
    return false;
}

bool mismatch_anchor(CReadAln* read_aln, char* md_string, int ref_start, bam1_t* bam_record) {
    if (!md_string) return false;
    
    // Create a copy of the MD string for parsing
    char* md_copy = Gstrdup(md_string);
    
    bool has_anchor_mismatch = (process_MD_string(read_aln, md_copy, ref_start) || 
                               process_cigar_string(read_aln, bam_record));
    
    GFREE(md_copy);
    return has_anchor_mismatch;
}

bool analyze_mismatches(CReadAln* read_aln, GSamRecord* brec, int bundle_start) {
    
    // Check left clipping proximity to splice site
    if (brec->clipL && 
        read_aln->segs[0].len() < junctionsupport + brec->clipL) {
        return true;
    }
    
    // Check right clipping proximity to splice site
    if (brec->clipR && 
        read_aln->segs.Last().len() < junctionsupport + brec->clipR) {
        return true;
    }
    
    // Check anchor mismatches
    if (mismatch_anchor(read_aln, brec->tag_str("MD"), bundle_start, brec->get_b())) {
        return true;
    }
    
    return false;
}

bool determine_mismatch_status(CReadAln* read_aln, GSamRecord* brec,
                           double nm, int bundle_start) {
    
    // Long reads assumed to have higher mismatch rate
    if (read_aln->longread) return true;
    
    // Check mismatch fraction threshold
    if (nm / read_aln->len > mismatchfrac) return true;
    
    // Detailed mismatch analysis for reads with any mismatches
    if (nm > 0) return analyze_mismatches(read_aln, brec, bundle_start);
    
    return false;
}

// a "well-anchored junction" refers to a splice junction that is supported 
// by a sufficient number of reads spanning across the exon-exon boundary, 
// making it a confident and reliable detection of alternative splicing
bool is_well_anchored(CReadAln* read_aln, int junction_index) {
    return (read_aln->segs[junction_index].len() > longintronanchor &&
            read_aln->segs[junction_index + 1].len() > longintronanchor);
}

void update_junction_statistics(CReadAln* read_aln, bool has_mismatch, 
                             int nh, float read_count) {
    
    for (int i = 0; i < read_aln->juncs.Count(); i++) {
        CJunction* junction = read_aln->juncs[i];
        
        // Update mismatch counter
        if (has_mismatch || nh > 2) junction->nm += read_count;
        
        // Update good match counter for well-anchored junctions
        if (is_well_anchored(read_aln, i)) junction->mm += read_count;
        
        // Update total read support
        junction->nreads += read_count;
    }
}

void process_junction_mismatches(CReadAln* read_aln, GSamRecord* brec,
                              double nm, float read_count, int bundle_start) {
    
    if (!read_aln->juncs.Count()) return;
    
    // Adjust mismatch count for paired reads and clipping
    double adjusted_nm = adjust_mismatch_count(brec, nm);
    
    // Determine if read has significant mismatches
    bool has_mismatch = determine_mismatch_status(read_aln, brec, adjusted_nm, bundle_start);
    
    // Update junction statistics
    update_junction_statistics(read_aln, has_mismatch, read_aln->nh, read_count);
}

std::string create_read_id(const char* read_name, int position, int hi) {
    std::string id(read_name);
    id += '-';
    id += position;
    id += ".=";
    id += hi;
    return id;
}

void add_pair_if_new(CReadAln* read_aln, int pair_index, float read_count) {
    // Check if pairing already exists
    for (int i = 0; i < read_aln->pair_idx.Count(); i++) {
        if (read_aln->pair_idx[i] == pair_index) {
            return; // Pairing already exists
        }
    }
    
    // Add new pairing
    read_aln->pair_idx.Add(pair_index);
    read_aln->pair_count.Add(read_count);
}

void establish_pairing(BundleData& bundle, int read1_index, int read2_index, 
                     float read_count) {
    
    // Add pairing information to both reads if not already present
    add_pair_if_new(bundle.reads[read1_index], read2_index, read_count);
    add_pair_if_new(bundle.reads[read2_index], read1_index, read_count);
}

void process_paired_reads(BundleData& bundle, int bundle_start, int read_start,
                       int read_index, GSamRecord* brec,
                       float read_count, int hi, GHash<int>& hashread) {
    
    // Only process pairs on same chromosome/contig
    if (brec->refId() != brec->mate_refId()) return;
    
    int mate_start = brec->mate_start();
    
    // Ignore pairs in previous bundles
    if (mate_start < bundle_start) return;
    
    // Handle pair processing based on mate position
    if (mate_start <= read_start) {
        std::string read_id = create_read_id(bundle.reads[read_index]->brec->name(),
                                      mate_start, hi);
    
    	const int* mate_index = hashread[read_id.c_str()];
    	if (mate_index) {
        	establish_pairing(bundle, read_index, *mate_index, read_count);
        	hashread.Remove(read_id.c_str());
    	}
        
    } else {
        std::string read_id = create_read_id(brec->name(), read_start, hi);
    	hashread.Add(read_id.c_str(), read_index);
    }
}

void process_read_in(uint bundle_start, uint bundle_end, BundleData& bundle,
                 GHash<int>& hashread, GReadAlnData& alndata, GSamRecord* brec) {
    
    // Skip secondary alignments
    if (brec->flags() & BAM_FSECONDARY) return;

    // TODO: collapse reads that have same start, end, and CIGAR
    // keep in list or something so attributes are maintained
    
    // Extract alignment information
    char strand = alndata.strand;
    int nh = alndata.nh;
    int hi = alndata.hi;
    int read_start = brec->start;
    double nm = static_cast<double>(brec->tag_int("NM"));
    float unitig_cov = brec->tag_float("YK");
    
    // Create new read alignment
    CReadAln* readaln = new CReadAln(strand, nh, brec->start, brec->end);
    readaln->longread = (longreads || brec->uval);
    
    // Process exons and junctions
    process_exons_and_junctions(brec, readaln, bundle.junction, alndata);
    
    // Add read to bundle
    int n = bundle.reads.Add(readaln);
    bundle.reads[n]->brec = brec;
    
    // Update bundle end if necessary
    update_bundle_end(bundle_end, bundle, brec->end);
    
    // Calculate read count with multi-mapping correction
    float read_count = calculate_read_count(brec, unitig_cov);
    bundle.reads[n]->read_count += read_count;
    
    // Process junction mismatches
    process_junction_mismatches(bundle.reads[n], brec, nm, read_count, bundle_start);
    
    // Handle paired-end reads
    process_paired_reads(bundle, bundle_start, read_start, n, brec, 
        read_count, hi, hashread);
}