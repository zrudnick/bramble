// reads.h

#ifndef READS_H
#define READS_H

#include "bramble.h"

bool has_deletion_at_junction(GSamRecord* brec, int junction_index);

bool has_guide_match(GReadAlnData& alndata, int junction_index);

CJunction* create_junction(GSamRecord* brec, int exon_index, 
        char strand, GReadAlnData& alndata, GList<CJunction>& junction);

void process_exons_and_junctions(GSamRecord* brec, CReadAln* readaln,
        GList<CJunction>& junction, GReadAlnData& alndata);

void update_bundle_end(uint bundle_end, BundleData& bundle, int read_end);

float calculate_read_count(GSamRecord* brec, float unitig_cov);

double adjust_mismatch_count(GSamRecord* brec, double nm);

bool process_matching_bases(char*& parser, uint& parsed_length);

void update_segment_tracking(CReadAln* read_aln, uint& segment_index, 
        uint& read_length, int parsed_length);

bool is_deletion_near_splice_site(CReadAln* read_aln, int segment_index, 
        int parsed_length, int deletion_length);

bool process_deletion(CReadAln* read_aln, char*& c, uint& segment_index, 
        uint& parsed_length, uint& read_length);

bool is_mismatch_near_splice_site(CReadAln* read_aln, int segment_index, 
        int parsed_length);

bool process_mismatch(CReadAln* read_aln, uint& segment_index, int parsed_length, 
        uint& read_length);

bool process_MD_string(CReadAln* read_aln, char* md_string, int ref_start);

bool is_insertion_near_splice_site(CReadAln* read_aln, int segment_index, int parsed_length);

bool process_insertion(CReadAln* read_aln, uint& segment_index, int parsed_length, 
        uint& read_length);

bool process_cigar_string(CReadAln* read_aln, bam1_t* bam_record);

bool mismatch_anchor(CReadAln* read_aln, char* md_string, int ref_start, 
        bam1_t* bam_record);

bool analyze_mismatches(CReadAln* read_aln, GSamRecord* brec, int bundle_start);

bool determine_mismatch_status(CReadAln* read_aln, GSamRecord* brec, 
        double nm, int bundle_start);

bool is_well_anchored(CReadAln* read_aln, int junction_index);

void update_junction_statistics(CReadAln* read_aln, bool has_mismatch, 
        int nh, float read_count);

void process_junction_mismatches(CReadAln* read_aln, GSamRecord* brec, double nm, 
        float read_count, int bundle_start);

std::string create_read_id(const char* read_name, int position, int hi);

void add_pair_if_new(CReadAln* read_aln, int pair_index, float read_count);

void establish_pairing(BundleData& bundle, int read1_index, int read2_index, 
        float read_count);

void process_paired_reads(BundleData& bundle, int bundle_start, int read_start,
        int read_index, GSamRecord* brec, float read_count, int hi, GHash<int>& hashread);

void process_read_in(uint bundle_start, uint bundle_end, BundleData& bundle,
        GHash<int>& hashread, GReadAlnData& alndata, GSamRecord* brec);

#endif