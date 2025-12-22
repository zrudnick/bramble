
#pragma once

#include "GArgs.h"
#include "GBitVec.h"
#include "GHashMap.hh"
#include "GSam.h"
#include "GStr.h"
#include "gff.h"
#include "time.h"

#include "types.h"

namespace bramble {

  struct BamIO;

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
    short int nh;
    uint32_t len;
    float read_count;             // keeps count for all reads (including paired and unpaired)
    bool unitig : 1;              // set if read come from an unitig
    bool longread : 1;            // set if read comes from long read data
    GVec<float> pair_count;       // keeps count for all paired reads
    GVec<int> pair_idx;           // keeps indices for all pairs in assembly mode, or all
                                    // reads that were collapsed in merge mode

    // For mate position calculations
    uint32_t mate_genomic_pos;
    tid_t mate_transcript_id;
    bool mate_found;

    GVec<GSeg> segs;  // "exons"
    GSamRecord *brec; // store BAM record

    CReadAln(char _strand = 0, short int _nh = 0, 
      int rstart = 0, int rend = 0)
      : GSeg(rstart, rend),
        strand(_strand), 
        nh(_nh), 
        len(0), 
        read_count(0), 
        unitig(false),
        longread(false), 
        pair_count(), 
        pair_idx(), 
        mate_genomic_pos(0),
        mate_transcript_id(), 
        mate_found(false), 
        segs(), 
        brec() {}

    CReadAln(CReadAln &rd) 
      : GSeg(rd.start, rd.end) { // copy contructor
          strand = rd.strand;
          nh = rd.nh;
          len = rd.len;
          read_count = rd.read_count;
          unitig = rd.unitig;
          longread = rd.longread;
          pair_count = rd.pair_count;
          pair_idx = rd.pair_idx;
        }

    ~CReadAln() { delete brec; }
  };

  // Bundle data structure, holds input data parsed from BAM file
  struct BundleData {
    BundleStatus status; // bundle status

    int idx; // index in the main bundles array
    int start; // bundle start
    int end; // bundle end
    unsigned long num_reads; // number of reads in this bundle

    GStr refseq;           // reference sequence name
    char* gseq;            // genomic sequence for the bundle
    GPVec<GffObj> guides;  // all guides in bundle
    std::vector<CReadAln *> reads; // all reads in bundle

    BundleData()
      : status(BundleStatus::BUNDLE_STATUS_CLEAR), 
        idx(0), 
        start(0), 
        end(0), 
        refseq(), 
        gseq(NULL), 
        guides(false),
        reads() {}
  
    // Called when the bundle is valid and ready to be processed
    void getReady(int curr_start, int curr_end) {
      start = curr_start;
      end = curr_end;
      status = BundleStatus::BUNDLE_STATUS_READY;
    }

    void keepGuide(GffObj *t) {
      guides.Add(t);
    }

    void Clear() {

      guides.Clear();

      // delete each CReadAln*
      for (auto &aln : reads) {
        delete aln;
      }
      reads.clear();

      // TODO: Think about managing this with a unique_ptr
      GFREE(gseq); 

      start = 0;
      end = 0;
      status = BundleStatus::BUNDLE_STATUS_CLEAR;
    }

    ~BundleData() { Clear(); }
  };

  // -------- function definitions

  void process_exons(GSamRecord *brec, CReadAln *readaln);

  int add_new_read(BundleData &bundle, CReadAln *readaln, GSamRecord *brec);

  void update_bundle_end(uint bundle_end, BundleData &bundle, int read_end);

  float calculate_read_count(GSamRecord *brec);

  std::string create_read_id(const char *read_name, int position, int hi);

  void add_pair_if_new(CReadAln *read_aln, int pair_index, float read_count);

  void establish_pairing(BundleData &bundle, int read1_index, int read2_index,
                      float read_count);

  void process_paired_reads(BundleData &bundle, int bundle_start, int read_start,
                          int read_index, GSamRecord *brec, float read_count,
                          int hi, GHash<int> &hashread);

  void process_read_in(uint bundle_start, uint bundle_end, BundleData &bundle,
                    GHash<int> &hashread, GSamRecord *brec, char strand, 
                    int nh, int hi);

  void build_bundles(BamIO* io, BundleData *bundles, BundleData* bundle, 
                    GPVec<BundleData> *bundle_queue, GVec<GRefData> refguides, 
                    int n_refguides, std::shared_ptr<GffNames> guide_seq_names);

}

