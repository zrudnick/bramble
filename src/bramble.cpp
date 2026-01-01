
#include <numeric>
#include <algorithm>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "CLI/CLI11.hpp"
#ifndef NOTHREADS
#include "GThreads.h"
#endif
#include "htslib/sam.h"
#include "quill/Backend.h"
#include "quill/Frontend.h"
#include "quill/LogMacros.h"
#include "quill/Logger.h"
#include "quill/sinks/ConsoleSink.h"

#include "types.h"
#include "bundles.h"
#include "bramble.h"
#include "evaluate.h"
#include "g2t.h"
#include "threads.h"

#include <IITree.h>
#include "GArgs.h"
#include "GBitVec.h"
#include "GSam.h"
#include "GStr.h"
#include "gff.h"
#include "time.h"

using namespace bramble;

#define VERSION "1.0.0"

#define USAGE                                                                  \
  "Bramble v" VERSION " usage:\n\n\
bramble <in.bam ..> [-G <guide_gff>] [-o <out.bam>] [-p <cpus>] [-S <genome.fa>] \n\
 [--help] [--version] [--verbose] [--long] [--fr] [--rf] \n\
 \n\
Project spliced genomic alignments into transcriptomic space.\n\
Options:\n\
 --help     : print this usage message and exit\n\
 --version  : print just the version at stdout and exit\n\
 --verbose  : verbose (log bundle processing details)\n\
 --long     : alignments are from long reads?\n\
 -G <file>  : reference annotation to use for guiding the BAM conversion (GTF/GFF)\n\
 -o <file>  : output path/file name for the projected alignments (default: stdout)\n\
 -p <int>   : number of threads (CPUs) to use (default: 1)\n\
 -S <file>  : genome sequence file (FASTA format)\n\
"

// --fr       : assume stranded library fw-firststrand\n\
// --rf       : assume stranded library fw-secondstrand\n\

bool VERBOSE = false;     // Verbose, --verbose
bool DEBUG_ = false;
bool LONG_READS = false;  // BAM file contains long reads, --long
bool FR_STRAND = true;   // Read 1 is on forward strand, --fr
bool RF_STRAND = false;   // Read 1 is on reverse strand, --fr
bool USE_FASTA = false;   // Use FASTA for reducing soft clips
bool SOFT_CLIPS = true;   // Add soft clips
bool STRICT = false;      // Use for strict boundary adherence

FILE *f_out = NULL;       // Default: stdout
uint8_t n_threads = 1;    // Threads, -p
GStr bam_path_out;        // Output BAM path, -o
GStr guide_gff;           // Reference annotation, -G
GFastaDb *gfasta = NULL;  // FASTA file, -S
const char *bam_file_in;  // Input BAM path
const char *bam_file_out; // Output BAM path
const char *out_dir;      // Output folder

#ifndef NOTHREADS
// Threading: single producer, multiple consumers
// Main thread/program is always loading the producer

GMutex data_mutex; // Manage availability of data records ready to be loaded by
                   // main thread
GVec<int> clear_data_pool; // Indices of data bundles cleared for loading by
                           // main thread (clear data pool)
GConditionVar
    have_bundles; // Will notify a thread that a bundle was loaded in the ready
                  // queue (or that no more bundles are coming)
char bundle_work =
    1; // Bit 0 set if bundles are still being prepared (BAM file not exhausted
       // yet) Bit 1 set if there are Bundles ready in the queue

GMutex wait_mutex;       // Controls threads_waiting (idle threads counter)
uint8_t threads_waiting; // Idle worker threads
GConditionVar
    have_threads; // Will notify the bundle loader when a thread
                  // Is available to process the currently loaded bundle

GConditionVar have_clear; // Will notify when bundle buf space available
GMutex queue_mutex;       // Controls bundle_queue and bundles access
GFastMutex log_mutex;     // Only when verbose - to avoid mangling the log output
GFastMutex reading_mutex;
GFastMutex bam_io_mutex;  // Protects BAM io
#endif

bool no_more_bundles = false;
uint32_t dropped_reads;
uint32_t unresolved_reads;

std::tuple<uint32_t, uint32_t> 
get_guide_coordinates(GffExon *g, char strand, uint32_t chr_end) {
  uint32_t g_start, g_end;
  // if (strand == '+') {
    g_start = g->start;
    g_end = (g->end + 1);
  // } else {
  //   g_start = chr_end - (g->end + 1);
  //   g_end = chr_end - g->start;  
  // }
  return std::make_tuple(g_start, g_end);
}

struct GuideStruct {
  GVec<GRefData> refguides;
  int n_refguides;
  std::shared_ptr<GffNames> guide_seq_names;
  GffReader *gffreader;
};

GuideStruct guide_stuff(FILE *header_file) {
  // Open GFF/GTF file
  FILE *f = fopen(guide_gff.chars(), "r");
  if (f == NULL)
    GError("Error: could not open reference annotation file (%s)!\n",
           guide_gff.chars());

  // GffReader: transcripts only, sort by location
  GffReader *gffreader = new GffReader(f, true, true); // loading only recognizable transcript features
  gffreader->setRefAlphaSorted(); // alphabetical sorting of RefSeq IDs
  gffreader->showWarnings(DEBUG_);

  // keep attributes, merge close exons, no exon attributes
  // merge_close_exons must be false for correct transcriptome header
  // construction!
  gffreader->readAll(true, false, false);

  int n_refguides = gffreader->gseqtable.Count();
  if (n_refguides == 0 || gffreader->gflst.Count() == 0) {
    GError("Error: could not any valid reference transcripts in %s (invalid "
           "GTF/GFF file?)\n",
           guide_gff.chars());
  }

  GVec<GRefData> refguides;        // vector with guides for each chromosome
  refguides.setCount(n_refguides); // maximum reference guide ID

  // Process each transcript for guide loading and header generation
  int last_refid = -1;
  for (int i = 0; i < gffreader->gflst.Count(); i++) {
    GffObj *guide = gffreader->gflst[i];

    // Chromosome switch
    if (last_refid != guide->gseq_id) {
      last_refid = guide->gseq_id;
    }

    // Sanity check: make sure there are no exonless "genes" or other
    if (guide->exons.Count() == 0) {
      if (DEBUG_)
        // LOG_WARNING(logger, "Warning: exonless GFF {} feature with ID {} found, added "
        //          "implicit exon {}-{}.\n",
        //          guide->getFeatureName(), guide->getID(), guide->start,
        //          guide->end);
      guide->addExon(guide->start, guide->end); // should never happen!
    }

    // Header generation: calculate spliced transcript length by summing exon
    // lengths
    int transcript_len = 0;
    for (int j = 0; j < guide->exons.Count(); j++) {
      if (guide->exons[j] != NULL) {
        transcript_len += guide->exons[j]->end - guide->exons[j]->start + 1;
      }
    }

    // Write to SAM header if we have valid transcript info
    auto id = guide->getID();
    if (id != NULL && transcript_len > 0) {
      fprintf(header_file, "@SQ\tSN:%s\tLN:%d\n", id, transcript_len);
    }

    GRefData &grefdata = refguides[guide->gseq_id];
    grefdata.add(gffreader, guide); // transcripts already sorted by location
  }

  auto guide_seq_names = std::shared_ptr<GffNames>(GffObj::names);

  fprintf(header_file, "@CO\tGenerated from GTF: %s\n", guide_gff.chars());
  fclose(header_file);

  return {refguides, n_refguides, guide_seq_names, gffreader};
}

std::shared_ptr<g2tTree> build_g2t_tree(GuideStruct guideStuff, BamIO *io) {
  auto g2t = std::make_shared<g2tTree>();

  GList<GffObj> *guides;
  for (uint32_t refid = 0; refid < guideStuff.refguides.Count(); refid++) {
    guides = &(guideStuff.refguides[refid].rnas);
    uint32_t ref_len = io->get_header()->target_len[refid];

    for (int j = 0; j < guides->Count(); j++) {
      GffObj *guide = (*guides)[j];
      char strand = guide->strand;
      const char *tid_string = guide->getID(); // we don't strictly need this
#ifndef NOTHREADS
    bam_io_mutex.lock();
#endif
      tid_t tid = g2t->insertTidString(tid_string, io);
#ifndef NOTHREADS
    bam_io_mutex.unlock();
#endif

      uint32_t cum_len = 0;
      int idx;
      int exon_count = guide->exons.Count();
      for (int k = 0; k < exon_count; k++) {

        if (strand == '-') {
          idx = exon_count - k - 1;
        } else idx = k;
        auto coords = get_guide_coordinates(guide->exons[idx], strand, ref_len);
        uint32_t g_start = std::get<0>(coords);
        uint32_t g_end = std::get<1>(coords);
        
        g2t->addInterval(refid, g_start, g_end, tid, (uint8_t)k, 
          cum_len, strand);
        cum_len += (g_end - g_start);
      }
    }

    g2t->addRefLen(refid, ref_len);
    g2t->indexTrees(refid);
  }

  return g2t;
}

void free_guides(GuideStruct guideStuff) {
  // Delete gffreader data
  for (int i = 0; i < guideStuff.gffreader->gflst.Count(); i++) {
    GffObj *guide = guideStuff.gffreader->gflst[i];
    delete guide;
  }
  delete guideStuff.gffreader;
}

char get_strand(GSamRecord* brec) {
  char strand = brec->spliceStrand(); // tagged strand gets priority

  // Set strand if stranded library
  if ((strand == '.') && (FR_STRAND || RF_STRAND)) {

    if (brec->isPaired()) {
      if (brec->pairOrder() == 1) { // first read in pair
        if ((RF_STRAND && brec->revStrand()) || (FR_STRAND && !brec->revStrand()))
          strand = '+';
        else strand = '-';
      } else { // second read in pair
        if ((RF_STRAND && brec->revStrand()) || (FR_STRAND && !brec->revStrand()))
          strand = '-';
        else strand = '+';
      }
    } else { // read isn't paired
      if ((RF_STRAND && brec->revStrand()) || (FR_STRAND && !brec->revStrand()))
        strand = '+';
      else strand = '-';
    }
  }
  return strand;
}

void process_reads(std::shared_ptr<g2tTree> g2t, BamIO *io, 
                  std::shared_ptr<GffNames> guide_seq_names) {
  std::shared_ptr<ReadEvaluator> evaluator;
  if (LONG_READS) {
    evaluator = std::make_shared<LongReadEvaluator>();
  } else {
    evaluator = std::make_shared<ShortReadEvaluator>();
  }

  dropped_reads = 0;
  unresolved_reads = 0;
  int all_reads = 0;

  bool more_alignments = true;
  uint32_t unmapped_reads = 0;

  GSamRecord *brec = nullptr;
  std::vector<CReadAln *> reads;
  GHash<int> hashread;

  read_id_t id = 0;
  uint32_t chunk_size = 50000;

  std::string prev_read_name;
  std::string read_name;
  bool new_read_name;

  auto flush = [&]() {
    convert_reads(reads, g2t, evaluator, io);
    for (auto &read : reads) {
      delete read;
    }
    reads.clear();
    id = 0;
  };

  while (more_alignments) { 
    if ((brec = io->next()) != NULL) {
      all_reads++;
      if (brec->isUnmapped()) {
        unmapped_reads++;
        continue;
      }

      read_name = brec->name();
      if (id == 0) {
        new_read_name = true;
      } else {
        new_read_name = (read_name != prev_read_name);
      }

      if (id > chunk_size && new_read_name) {
        flush();
      }

      id += 1;
      char strand = get_strand(brec);
      const char *ref_name = brec->refName();
      refid_t refid = guide_seq_names->gseqs.addName(ref_name);
      process_read_in(reads, id, brec, strand, refid, 
        brec->tag_int("NH"), brec->tag_int("HI"), hashread);

    } else {
      // No more alignments
      more_alignments = false;
      flush();
    }

    prev_read_name = read_name;
  }

  GMessage("#total input reads is:    %d\n", all_reads);
  GMessage("#unmapped reads is: %d\n", unmapped_reads);
  GMessage("#dropped reads is: %d\n", dropped_reads);
  GMessage("#unresolved reads is: %d\n\n", unresolved_reads);
}

int main(int argc, char *argv[]) {

  //quill::Backend::start();
  // Get a logger
  // This segfaults on my system, not sure why
  // quill::Logger *logger = quill::Frontend::create_or_get_logger(
  //     "root",
  //     quill::Frontend::create_or_get_sink<quill::ConsoleSink>("sink_id_1"));

  CLI::App app{
      "Project spliced genomic alignments into transcriptomic space",
      "Bramble"};

  std::string gff;
  std::string bam_file;
  std::string input_bam;
  app.add_option("in.bam", input_bam, "input bam file")->required();
  // app.add_flag("--fr", FR_STRAND, "assume stranded library fw-firststrand");
  // app.add_flag("--rf", RF_STRAND, "assume stranded library fw-secondstrand");
  app.add_flag("--verbose", VERBOSE, "verbose (log processing details)");
  app.add_flag("--long", LONG_READS, "alignments are from long reads?");
  app.add_option("-G", gff,
                 "reference annotation to use for guiding the assembly process "
                 "(GTF/GFF)")
      ->check(CLI::ExistingFile);
  app.add_option("-o", bam_file,
                 "output path/file name for the projected alignments")
      ->required();
  app.add_option("-p", n_threads, "number of threads (CPUs) to use")
      ->default_val(1);
  auto version_callback = [](int count) {
    std::cout << "version: " << VERSION << "\n";
    std::exit(EXIT_SUCCESS);
  };
  app.add_flag_function("-V,--version", version_callback, "print version.");

  CLI11_PARSE(app, argc, argv);

  // LOG_INFO(logger, "command line parsed --- running bramble");
  bam_file_in = input_bam.c_str();
  guide_gff = gff.c_str();
  bam_file_out = bam_file.c_str();

  // Create SAM header file
  const GStr header_path = "tx_header.sam";
  FILE *header_file = fopen(header_path.chars(), "w");
  if (header_file == NULL)
    GError("Error creating file: %s\n", header_path.chars());
  fprintf(header_file, "@HD\tVN:1.0\tSO:coordinate\n");

  if (VERBOSE || DEBUG_) {
    // print logo
    GMessage("   __                  __   __   \n");
    GMessage("  / /  _______ ___ _  / /  / /__ \n");
    GMessage(" / _ \\/ __/ _ `/  ' \\/ _ \\/ / -_)\n");
    GMessage("/_.__/_/  \\_,_/_/_/_/_.__/_/\\__/ \n\n");
    GMessage("## Running Bramble version %s ##\n\n", VERSION);
    GMessage("-- Loading reference annotation...\n");
  }
    
  GuideStruct guideStuff = guide_stuff(header_file);

  if (VERBOSE || DEBUG_) {
    GMessage("Reference annotation loaded! We found %d unique transcripts\n\n", guideStuff.gffreader->gflst.Count());
  
    if (LONG_READS) {
      GMessage("-- Using LONG READS mode\n");
    }
    else {
      GMessage("-- Using SHORT READS mode\n");
      GMessage("Have long-read input data? Try running with the --long flag\n"); 
    }
  }
  
  BamIO *io = new BamIO(bam_file_in, bam_file_out, header_path);
  io->start();

  auto g2t = build_g2t_tree(guideStuff, io);
  free_guides(guideStuff);
  process_reads(g2t, io, guideStuff.guide_seq_names);

  io->stop(); // close BAM reader & writer
  delete io;
  delete gfasta;
  
  f_out = stdout;
  fprintf(f_out, "# Bramble version %s\n", VERSION);
}
