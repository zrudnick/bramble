// reads.h

#ifndef READS_H
#define READS_H

#include "bramble.h"

int add_new_read(BundleData& bundle, CReadAln* readaln, GSamRecord* brec);

void process_exons(GSamRecord* brec, CReadAln* readaln, GReadAlnData& alndata);

void update_bundle_end(uint bundle_end, BundleData& bundle, int read_end);

float calculate_read_count(GSamRecord* brec);

std::string create_read_id(const char* read_name, int position, int hi);

void add_pair_if_new(CReadAln* read_aln, int pair_index, float read_count);

void establish_pairing(BundleData& bundle, int read1_index, int read2_index, 
        float read_count);

void process_paired_reads(BundleData& bundle, int bundle_start, int read_start,
        int read_index, GSamRecord* brec, float read_count, int hi, GHash<int>& hashread);

void process_read_in(uint bundle_start, uint bundle_end, BundleData& bundle,
        GHash<int>& hashread, GReadAlnData& alndata, GSamRecord* brec);

#endif

