#include <numeric>
#include <algorithm>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#ifdef USE_OPENMP
#include <omp.h>
#endif

#ifndef NOTHREADS
#include "GThreads.h"
#endif

#include "CLI/CLI11.hpp"
#include "htslib/sam.h"
#include "quill/Backend.h"
#include "quill/Frontend.h"
#include "quill/LogMacros.h"
#include "quill/Logger.h"
#include "quill/sinks/ConsoleSink.h"

#include "types.h"
#include "bundles.h"
#include "bramble.h"
#include "g2t.h"

#include "GArgs.h"
#include "GBitVec.h"
#include "GSam.h"
#include "GStr.h"
#include "gff.h"
#include "time.h"

#include "IITree.h"
#include <assert.h>

using namespace bramble;

// Global variables needed by bundling code
bool VERBOSE = false;
bool BRAMBLE_DEBUG = false;
bool LONG_READS = false;
bool FR_STRAND = false;
bool RF_STRAND = false;
bool USE_FASTA = false;
bool SOFT_CLIPS = false;
bool STRICT = false;

uint8_t n_threads = 1;
GFastaDb *gfasta = nullptr;


bool no_more_bundles = false;

// pseudocode
// load gtf/gff into GffReader
// sort
// get number of references - initiate array of IITrees, one per reference
// create guide bundles - iterate over sorted gffreader and create bundles of overlapping transcripts
// bundles are constructed in parallel with a lock on insert into the IITree
// iterate over reads - look up g2t tree in IITrees and process as before

struct TX {
    public:
        TX() = default;
        ~TX() = default;
        TX(int seqid, GffObj* tx) {
            this->seqid = seqid;
            this->tid = tx->getID();
            this->strand = tx->strand;
            this->start = tx->start;
            this->end = tx->end;

            // add exons
            for (int i = 0; i < tx->exons.Count(); i++) {
                exons.push_back(std::make_pair(tx->exons[i]->start, tx->exons[i]->end));
            }
        };

        int seqid;
        char strand;
        int start;
        int end;
        std::string tid;
        std::vector<std::pair<uint, uint>> exons;
};

struct SimpleBundle {
    int32_t start;
    int32_t end;
    std::vector<TX> guides;

    SimpleBundle() : start(-1), end(-1) {}

    void add_guide(const TX& tx) {
        if (guides.empty()) {
            start = tx.start;
            end = tx.end;
        } else {
            if (tx.start < start) start = tx.start;
            if (tx.end > end) end = tx.end;
        }
        guides.push_back(tx);
    }

    bool overlaps(const TX& tx) const {
        if (guides.empty()) return false;
        // Check if guide overlaps with current bundle range (with RUNOFF_DIST)
        return tx.start <= end + RUNOFF_DIST;
    }

    void clear() {
        start = -1;
        end = -1;
        guides.clear();
    }

    bool empty() const {
        return guides.empty();
    }
};

std::unique_ptr<g2tTree> make_g2t_tree(SimpleBundle& bundle, BamIO* io) {
    if (bundle.guides.size() == 0)
      return nullptr;

    auto g2t = std::make_unique<g2tTree>();

    // Insert all guide exons into the interval tree
    for (int i = 0; i < bundle.guides.size(); i++) {
        TX& guide = bundle.guides[i];

        std::string tid_string = guide.tid;
        const char* tid_cstr = tid_string.c_str();

        tid_t tid = g2t->insertTidString(tid_cstr, io);

        for (auto exon : guide.exons) {
            uint g_start, g_end;
            if (guide.strand == '+') {
                g_start = exon.first - bundle.start;
                g_end = (exon.second + 1) - bundle.start;
            }
            else {
                g_start = bundle.end - (exon.second + 1);
                g_end = bundle.end - exon.first;
            }

            g2t->addGuideExon(g_start, g_end, tid, guide.strand);
      }
    }

    g2t->buildAllTidChains();
    g2t->precomputeAllCumulativeLengths();

    return g2t;
}

int main(int argc, char *argv[])
{
    assert (argc == 4);

    // get gtf from args
    FILE *f = fopen(argv[1], "r");
    if (f == NULL){
        std::cerr << "Error: could not open reference annotation file (" << argv[1] << ")!" << std::endl;
        return -1;
    }

    FILE *f2 = fopen(argv[2], "r");
    if (f2 == NULL){
        std::cerr << "Error: could not open input bam file (" << argv[2] << ")!" << std::endl;
        return -1;
    };
    fclose(f2);

    // output bam init
    const char *bam_file_out = "out.bam";
    const GStr header_path = "tx_header.sam";

    // get bam file from second argument
    BamIO *io = new BamIO(argv[2], bam_file_out, header_path);
    io->start();

    // read all references from the input bam header
    bam_hdr_t *bam_header = io->get_header();
    int n_seqids = bam_header->n_targets;
    std::unordered_map<std::string, int> seqid2id;
    for (int i = 0; i < n_seqids; i++) {
        std::string ref_name = std::string(bam_header->target_name[i]);
        seqid2id[ref_name] = i;
    }

    // GffReader: transcripts only, sort by location
    GffReader *gffreader = new GffReader(f, true, true); // loading only recognizable transcript features
    gffreader->setRefAlphaSorted(); // alphabetical sorting of RefSeq IDs
    gffreader->readAll(true, false, false);

    int gtf_n_seqids = gffreader->gseqtable.Count();
    if (gtf_n_seqids == 0 || gffreader->gflst.Count() == 0) {
        std::cerr << "Error: could not find any valid reference sequences" << argv[1] << " (invalid GTF/GFF file?)" << std::endl;
        return -1;
    }

    // create an array to hold IITrees for each reference sequence
    std::vector<IITree<int32_t, g2tTree*>> iitrees(gtf_n_seqids);

    SimpleBundle current_bundle;
    int last_tid = -1;
    for (int i = 0; i < gffreader->gflst.Count(); i++) {
        GffObj *pGffObj = gffreader->gflst[i];
        std::string tx_seqid = pGffObj->getGSeqName();
        // get the numerical bam id for this reference sequence - fail if not found
        auto it = seqid2id.find(tx_seqid);
        if (it == seqid2id.end()) {
            std::cerr << "Error: could not find reference sequence " << tx_seqid << " from GTF/GFF in BAM file!" << std::endl;
            return -1;
        }
        int bam_id = it->second;

        TX tx(bam_id, pGffObj);

        // if this is a new reference sequence, or no overlap - flush the current bundle
        if (bam_id != last_tid || !current_bundle.overlaps(tx)) {
            // add current bundle to IITree
            if (!current_bundle.empty()) {
                iitrees[last_tid].add(current_bundle.start, current_bundle.end, make_g2t_tree(current_bundle, io).release());
                current_bundle.clear();
            }
        }
        current_bundle.add_guide(tx);
        last_tid = bam_id;
    }

    // flush any remaining bundle
    if (!current_bundle.empty()) {
        iitrees[last_tid].add(current_bundle.start, current_bundle.end, make_g2t_tree(current_bundle, io).release());
        current_bundle.clear();
    }

    // idex all trees
    for (size_t ref_id = 0; ref_id < iitrees.size(); ref_id++) {
        iitrees[ref_id].index();
    }

    GSamRecord *brec = nullptr;
    bool more_alignments = true;

    while (more_alignments) {
        // Get next alignment
        if ((brec = io->next()) != NULL) {
            if (brec->isUnmapped()) {
                continue;
            }
            std::cout<< "Read: " << brec->name() << " (" << brec->start << ", " << brec->end << ")" << std::endl;
        }
        else{
            more_alignments = false;
            continue;
        }
    }
    
    // std::vector<size_t> results;
    
    // tree.overlap(150, 160, results);

    // for (size_t idx : results) {
    //     std::cout << "    - [" << tree.start(idx) << ", " << tree.end(idx) << ") -> " 
    //               << tree.data(idx) << std::endl;
    // }

    io->stop(); // close BAM reader & writer
    delete io;

    // Clean up g2tTree pointers
    for (size_t ref_id = 0; ref_id < iitrees.size(); ref_id++) {
        for (size_t i = 0; i < iitrees[ref_id].size(); i++) {
            delete iitrees[ref_id].data(i);
        }
    }
    
    return 0;
}
