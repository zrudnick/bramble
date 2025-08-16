
// tree.h

#ifndef TREE_H_
#define TREE_H_

#include <algorithm>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "bramble.h"

void print_tree(g2tTree* g2t);

std::unique_ptr<g2tTree> make_g2t_tree(BundleData *bundle, BamIO *io);

std::string reverse_complement(const std::string &seq);

std::string extract_sequence(const char *gseq, uint start, uint length,
                             char strand);

bool check_backward_overhang(IntervalNode *interval, uint exon_start,
                             char read_strand, std::set<std::string> &exon_tids,
                             g2tTree *g2t, BundleData *bundle);

std::set<tid_t> 
collapse_intervals(std::vector<IntervalNode *> sorted_intervals,
                   uint exon_start, bool is_first_exon, char read_strand,
                   g2tTree *g2t, BundleData *bundle,
                   bool &used_backwards_overhang, uint &soft_clip,
                   IntervalNode *prev_last_interval);

bool check_forward_overhang(IntervalNode *interval, uint exon_end,
                            char read_strand, std::set<tid_t> &exon_tids,
                            g2tTree *g2t, BundleData *bundle);

uint get_match_pos(IntervalNode *interval, tid_t tid, char read_strand,
                   g2tTree *g2t, uint exon_start, bool used_backwards_overhang);

void process_read_out(BundleData *&bundle, 
                      uint read_index, 
                      g2tTree *g2t,
                      std::map<read_id_t, ReadInfo *> &bam_info,
                      std::vector<uint32_t> group);

void add_mate_info(const std::set<tid_t> &final_transcripts,
                   const std::set<tid_t> &read_transcripts,
                   const std::set<tid_t> &mate_transcripts,
                   const std::map<tid_t, uint> &read_positions,
                   const std::map<tid_t, uint> &mate_positions,
                   std::map<read_id_t, ReadInfo *> &read_info, 
                   std::map<bam_id_t, BamInfo *> &bam_info, uint read_index,
                   uint mate_index, uint mate_case);

void update_read_matches(ReadInfo *read_info,
                         const std::set<tid_t> &final_transcripts);

void process_mate_pairs(BundleData *bundle, 
                        std::map<read_id_t, ReadInfo *> &read_info,
                        std::map<bam_id_t, BamInfo *> &bam_info);

void free_read_data(AlnGroups* aln_groups, std::map<read_id_t, ReadInfo *> &read_info,
                    std::map<bam_id_t, BamInfo *> &bam_info);

void convert_reads(BundleData *bundle, BamIO *io);

#endif
