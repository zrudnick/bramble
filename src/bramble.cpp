// Bramble v1.0.0

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
#include "threads.h"

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
bool DEBUG = false;
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
  auto input_opt =
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

  if (VERBOSE || DEBUG) {
    // print logo
    GMessage("   __                  __   __   \n");
    GMessage("  / /  _______ ___ _  / /  / /__ \n");
    GMessage(" / _ \\/ __/ _ `/  ' \\/ _ \\/ / -_)\n");
    GMessage("/_.__/_/  \\_,_/_/_/_/_.__/_/\\__/ \n\n");
    GMessage("## Running Bramble version %s ##\n\n", VERSION);
    GMessage("-- Loading reference annotation...\n");
  }
    

  // Open GFF/GTF file
  FILE *f = fopen(guide_gff.chars(), "r");
  if (f == NULL)
    GError("Error: could not open reference annotation file (%s)!\n",
           guide_gff.chars());

  // GffReader: transcripts only, sort by location
  GffReader *gffreader = new GffReader(f, true, true); // loading only recognizable transcript features
  gffreader->setRefAlphaSorted(); // alphabetical sorting of RefSeq IDs
  gffreader->showWarnings(DEBUG);

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
      if (DEBUG)
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

  if (VERBOSE || DEBUG) {
    GMessage("Reference annotation loaded! We found %d unique transcripts\n\n", gffreader->gflst.Count());
  
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

#ifndef NOTHREADS

#define DEF_TSTACK_SIZE 8388608
  size_t def_stack_size = DEF_TSTACK_SIZE;

#ifdef _GTHREADS_POSIX_
  size_t tstack_size = GThread::defaultStackSize();
  if (tstack_size < DEF_TSTACK_SIZE)
    def_stack_size = DEF_TSTACK_SIZE;
  if (DEBUG) {
    if (tstack_size < def_stack_size) {
      GMessage("Default stack size for threads: %d (increased to %d)\n",
               tstack_size, def_stack_size);
    } else
      GMessage("Default stack size for threads: %d\n", tstack_size);
  }
#endif

  GThread *threads = new GThread[n_threads]; // threads for processing bundles
  GPVec<BundleData> bundle_queue(false);
  BundleData *bundles =
      new BundleData[n_threads + 1];         // redef with more bundles

  clear_data_pool.setCapacity(n_threads + 1);

  WorkerArgs worker_args(&bundle_queue, io);

  // Start worker threads
  for (int b = 0; b < n_threads; b++) {
    threads[b].kickStart(worker_thread, (void *)&worker_args, def_stack_size);
    bundles[b + 1].idx = b + 1;
    clear_data_pool.Push(b);
  }
  BundleData* bundle = &(bundles[n_threads]);

#else

  // Just put everything into the same bundle
  BundleData bundles[1];
  BundleData *bundle = &(bundles[0]);
#endif

  // Build bundles and call core functions
  build_bundles(io, bundles, bundle, &bundle_queue, refguides, 
    n_refguides, guide_seq_names);

  // Delete thread and bundle arrays
  #ifndef NOTHREADS
    for (int t = 0; t < n_threads; t++) {
      threads[t].join();
    }
    if (DEBUG) {
      // LOG_INFO(logger, " All threads finished.\n");
    }
    delete[] threads;
    delete[] bundles;
  #endif

  io->stop(); // close BAM reader & writer
  delete io;
  delete gfasta;

  // Delete gffreader data
  for (int i = 0; i < gffreader->gflst.Count(); i++) {
    GffObj *guide = gffreader->gflst[i];
    delete guide;
  }
  delete gffreader;
  
  f_out = stdout;
  fprintf(f_out, "# Bramble version %s\n", VERSION);
}
