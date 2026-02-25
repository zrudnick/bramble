
#pragma once
#include <algorithm>
#include <iostream>
#include <memory>
#include <vector>
#include <functional>

#include "GArgs.h"
#include "GBitVec.h"
#include "GHashMap.hh"
#include "GSam.h"
#include "GStr.h"
#include "gff.h"
#include "time.h"
#include "types.h"

#include "quill/Backend.h"
#include "quill/Frontend.h"
#include "quill/LogMacros.h"
#include "quill/Logger.h"
#include "quill/sinks/ConsoleSink.h"

namespace bramble {

  struct g2tTree;
  struct ReadEvaluator;

  struct BamIO {
  protected:
    std::unique_ptr<GSamReader> reader{nullptr};
    std::unique_ptr<GSamWriter> writer{nullptr};
    std::shared_ptr<GSamRecord> current_record{nullptr};
    GStr input_bam;
    GStr output_bam;
    GStr header_sam;

  public:
    BamIO(const GStr &in_bam, 
          const GStr &out_bam, 
          const GStr &sam_header)
      : input_bam(in_bam), 
        output_bam(out_bam), 
        header_sam(sam_header) {}

    void start() {
      reader = std::make_unique<GSamReader>(input_bam.chars());
      writer = std::make_unique<GSamWriter>(output_bam.chars(), 
        header_sam.chars());

      if (GSamRecord* next_rec = reader->next()) { 
        current_record = std::shared_ptr<GSamRecord>(
          next_rec, [](GSamRecord* p) { delete p; }); 
      }
    }

    std::shared_ptr<GSamRecord>next() {
      auto result = current_record; 
      if (GSamRecord* next_rec = reader->next()) { 
        current_record = std::shared_ptr<GSamRecord>(
          next_rec, [](GSamRecord* p){ delete p; }); 
      } else { 
        current_record.reset(); 
      } 
      return result;
    }

    void write(bam1_t *b) {
      if (writer && b != nullptr) {
        writer->write(b);
      }
    }

    int32_t get_tid(const char *transcript_id) {
      return writer->get_tid(transcript_id);
    }

    bam_hdr_t* get_header() { 
      return reader->header(); 
    }

    void stop() { 
      current_record.reset(); 
    }
  };

  // Collect all refguide transcripts for a single genomic sequence
  struct GRefData {
    GList<GffObj> rnas; // all transcripts on this genomic seq
    int gseq_id;
    const char *gseq_name;
    GRefData(int gid = -1)
        : rnas(false, false, false), gseq_id(gid), gseq_name(NULL) {
      gseq_id = gid;
      if (gseq_id >= 0)
        gseq_name = GffObj::names->gseqs.getName(gseq_id);
    }

    void add(GffReader *gffr, GffObj *t) {
      if (gseq_id >= 0) {
        if (gseq_id != t->gseq_id)
          GError("Error: invalid call to GRefData::add() - different genomic "
                "sequence!\n");
      } else { // adding first transcript, initialize storage
        gseq_id = t->gseq_id;
        gseq_name = t->getGSeqName();
        if (gffr->gseqtable[gseq_id] == NULL)
          GError("Error: invalid genomic sequence data (%s)!\n", gseq_name);
        rnas.setCapacity(gffr->gseqtable[gseq_id]->fcount);
      }
      rnas.Add(t);
      t->isUsed(true);
      // setLocus(t); //use the GRefLocus::mexons to quickly find an overlap with
      // existing loci, or create a new one
    }

    bool operator==(GRefData &d) { return gseq_id == d.gseq_id; }
    bool operator<(GRefData &d) { return (gseq_id < d.gseq_id); }
  };

  // Bundle status
  enum BundleStatus {
    BUNDLE_STATUS_CLEAR = 0,      // Available for loading/prepping
    BUNDLE_STATUS_LOADING = 1,    // Being prepared by the main thread (there can be only one)
    BUNDLE_STATUS_READY = 2       // Ready to be processed, or being processed
  };

  // Read Alignment
  struct CReadAln : public GSeg {
    char strand;                  // 1, 0 (unknown), -1 (reverse)
    refid_t refid;                // chromosome origin
    float read_count;             // keeps count for all reads (including paired and unpaired)
    GVec<float> pair_count;       // keeps count for all paired reads
    GVec<int> pair_idx;           // keeps indices for all paired reads
    GVec<GSeg> segs;  // "exons"
    std::shared_ptr<GSamRecord> brec; // store BAM record

    CReadAln(char _strand = 0, refid_t id = 0, 
      int rstart = 0, int rend = 0)
      : GSeg(rstart, rend),
        strand(_strand), 
        refid(id),
        read_count(0),
        pair_count(), 
        pair_idx(), 
        segs(), 
        brec() {}

    ~CReadAln() {}
  };

  // Bundle data structure, holds input data parsed from BAM file
  struct BundleData {
    BundleStatus status; // bundle status
    int idx; // index in the main bundles array
    std::vector<CReadAln> reads; // all reads in bundle

    quill::Logger *logger;
    std::shared_ptr<g2tTree> g2t; 
    std::shared_ptr<ReadEvaluator> evaluator;

    BundleData()
      : status(BundleStatus::BUNDLE_STATUS_CLEAR), 
        idx(0),
        reads(),
        logger(),
        g2t(),
        evaluator() {}
  
    // Called when the bundle is valid and ready to be processed
    void getReady() {
      status = BundleStatus::BUNDLE_STATUS_READY;
    }

    void Clear() {
      reads.clear();
      status = BundleStatus::BUNDLE_STATUS_CLEAR;
    }

    ~BundleData() { 
      Clear(); 
    }
  };

  struct WorkerArgs {
    GPVec<BundleData> *bundle_queue;
    BamIO *io;

    WorkerArgs(GPVec<BundleData> *q, BamIO *d) 
      : bundle_queue(q), 
        io(d) {}
  };

}
