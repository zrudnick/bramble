
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

#include <IITree.h>
#include "GArgs.h"
#include "GBitVec.h"
#include "GSam.h"
#include "GStr.h"
#include "gff.h"
#include "time.h"

#include "types.h"
#include "bam.h"
#include "bramble.h"
#include "evaluate.h"
#include "g2t.h"
#include "threads.h"
#include "sw.h"

using namespace bramble;

#define VERSION "1.0.0"

#define USAGE                                                                  \
  "Bramble v" VERSION " usage:\n\n\
bramble <in.bam ..> [-G <guide_gff>] [-o <out.bam>] [-p <cpus>] [-S <genome.fa>] \n\
 [--help] [--version] [--quiet] [--long] [--fr] [--rf] \n\
 \n\
Project spliced genomic alignments into transcriptomic space.\n\
Options:\n\
 --help     : print this usage message and exit\n\
 --version  : print just the version at stdout and exit\n\
 --quiet    : turn off verbose mode (log bundle processing details)\n\
 --long     : alignments are from long reads?\n\
 -G <file>  : reference annotation to use for guiding the BAM conversion (GTF/GFF)\n\
 -o <file>  : output path/file name for the projected alignments (default: stdout)\n\
 -p <int>   : number of threads (CPUs) to use (default: 1)\n\
 -S <file>  : genome sequence file (FASTA format)\n\
"

// --fr       : assume stranded library fw-firststrand\n\
// --rf       : assume stranded library fw-secondstrand\n\

bool QUIET = false;     // Turn off verbose mode, --quiet
bool BRAMBLE_DEBUG = false;
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
uint32_t total_reads;
uint32_t unmapped_reads;
uint32_t dropped_reads;
uint32_t unresolved_reads;

std::shared_ptr<g2tTree> build_g2t_tree(GVec<GRefData> refguides, 
                                        int n_refguides, BamIO *io) {
  auto g2t = std::make_shared<g2tTree>();

  GList<GffObj> *guides;

  std::vector<IntervalData> intervals;
  intervals.reserve(100);
  
  for (uint32_t refid = 0; refid < refguides.Count(); refid++) {
    guides = &(refguides[refid].rnas);
    const char* ref_name = refguides[refid].gseq_name;

    for (int j = 0; j < guides->Count(); j++) {
      GffObj *guide = (*guides)[j];
      char strand = guide->strand;
      const char *tid_string = guide->getID();
      // todo: validate that these are contiguous
      
#ifndef NOTHREADS
      bam_io_mutex.lock();
#endif
      tid_t tid = g2t->insertTidString(tid_string, io);
#ifndef NOTHREADS
      bam_io_mutex.unlock();
#endif

      std::string tid_string_print = g2t->getTidName(tid);

      int exon_count = guide->exons.Count();
      intervals.clear();
      uint32_t pos_start = 0;

      for (int k = 0; k < exon_count; k++) {
        int idx = (strand == '-') ? (exon_count - k - 1) : k;
        
        uint32_t g_start = guide->exons[idx]->start;
        uint32_t g_end = guide->exons[idx]->end + 1;

        IntervalData interval;
        interval.start = g_start;
        interval.end = g_end;
        interval.idx = (uint8_t)idx;
        interval.pos_start = pos_start;
        intervals.push_back(interval);
        
        pos_start += g_end - g_start;
      }

      uint32_t transcript_len = pos_start;
      for (int k = 0; k < exon_count; k++) {
        auto interval = intervals[k];

        if (k > 0) {
          interval.prev_start = intervals[k-1].start;
          interval.prev_end = intervals[k-1].end;
          interval.has_prev = true;
        } else {
          interval.prev_start = 0;
          interval.prev_end = 0;
          interval.has_prev = false;
        }

        if (k < exon_count - 1) {
          interval.next_start = intervals[k+1].start;
          interval.next_end = intervals[k+1].end;
          interval.has_next = true;
        } else {
          interval.next_start = 0;
          interval.next_end = 0;
          interval.has_next = false;
        }

        interval.transcript_len = transcript_len;
        g2t->addInterval(refid, tid, interval, strand, 
          ref_name);
      }
    }

    g2t->indexTrees(refid);
  }

  return g2t;
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

void process_exons(GSamRecord *brec, CReadAln *read, const char *ref_name) {
  auto exons = brec->exons;
  for (int i = 0; i < exons.Count(); i++) {
    auto exon = exons[i];
    exon.end++; // switch to exclusive end

    read->len += exon.len();
    read->segs.Add(exon);
  }
}

float calculate_read_count(GSamRecord *brec) {
  float unitig_cov = brec->tag_float("YK");
  float read_count = static_cast<float>(brec->tag_int("YC"));
  if (read_count == 0)
    read_count = 1;
  if (unitig_cov)
    read_count = unitig_cov;
  return read_count;
}

std::string 
create_read_id(const char *read_name, int32_t pos, int32_t hi) {
  std::string id(read_name);
  id += '-'; id += pos; id += ".="; id += hi;
  return id;
}

void add_pair_if_new(CReadAln *read, int pair_id, float read_count) {
  for (int i = 0; i < read->pair_idx.Count(); i++) {
    if (read->pair_idx[i] == pair_id) {
      return; // pairing already exists
    }
  }

  read->pair_idx.Add(pair_id);
  read->pair_count.Add(read_count);
}

void process_pairs(std::vector<CReadAln *> &reads, 
                  read_id_t id, GSamRecord *brec, int32_t hi, 
                  float read_count, GHash<int> &hashread) {

  // Only process pairs on same chromosome/contig
  if (brec->refId() != brec->mate_refId())
    return;

  int32_t read_start = reads[id]->start;
  int32_t mate_start = brec->mate_start();

  if (mate_start <= read_start) {
    std::string read_id =
      create_read_id(reads[id]->brec->name(), mate_start, hi);
    int *mate_id = hashread[read_id.c_str()];
    if (mate_id) {
      add_pair_if_new(reads[id], *mate_id, read_count);
      add_pair_if_new(reads[*mate_id], id, read_count);
      hashread.Remove(read_id.c_str());
    }
  } else {
    std::string read_id = create_read_id(brec->name(), read_start, hi);
    hashread.Add(read_id.c_str(), id);
  }
}

void process_read_in(std::vector<CReadAln *> &reads, read_id_t id,
                    GSamRecord *brec, char strand, refid_t refid, 
                    const char *ref_name, int32_t nh, int32_t hi, GHash<int> &hashread) {

  // Create new read alignment
  CReadAln *read = new CReadAln(strand, refid, nh, brec->start, brec->end);
  read->brec = brec;
  reads.emplace_back(read);
  
  process_exons(brec, read, ref_name);
  float read_count = calculate_read_count(brec);
  process_pairs(reads, id, brec, hi, read_count, hashread);
}

void process_reads(std::shared_ptr<g2tTree> g2t, BamIO *io,  
                  BundleData *bundles, BundleData* bundle, 
                  GPVec<BundleData> *bundle_queue) {

  std::shared_ptr<ReadEvaluator> evaluator;
  if (LONG_READS) {
    evaluator = std::make_shared<LongReadEvaluator>();
  } else {
    evaluator = std::make_shared<ShortReadEvaluator>();
  }

  auto guide_seq_names = std::shared_ptr<GffNames>(GffObj::names);

  unmapped_reads = 0;
  total_reads = 0;
  dropped_reads = 0;
  unresolved_reads = 0;

  bool more_alignments = true;

  GSamRecord *brec = nullptr;
  std::vector<CReadAln *> reads;
  GHash<int> hashread;

  read_id_t id = 0;
  uint32_t chunk_size = 1000000;

  std::string prev_read_name;
  std::string read_name;
  bool new_read_name;

  while (more_alignments) { 
    bool new_bundle = false;
    char strand;
    const char *ref_name;
    refid_t refid;

    if ((brec = io->next()) != NULL) {
      total_reads++;
      if (brec->isUnmapped()) {
        unmapped_reads++;
        continue;
      }

      read_name = brec->name();
      strand = get_strand(brec);
      ref_name = brec->refName();
      refid = guide_seq_names->gseqs.addName(ref_name);

      if (id == 0) {
        new_read_name = true;
      } else {
        if (read_name != prev_read_name) {
          new_read_name = true;
        } else {
          new_read_name = false;
        }
      }

      if (id >= chunk_size && new_read_name) {
        new_bundle = true;
      }

    } else {
      // No more alignments
      more_alignments = false;
      new_bundle = true;
    }

    if (new_bundle) {
      hashread.Clear();
      bundle->g2t = g2t;
      bundle->evaluator = evaluator;

      // Push completed bundle
      push_bundle(bundle, bundle_queue);
      id = 0;

      if (!more_alignments) {
        wait_for_bundles();
        break;
      }

#ifndef NOTHREADS
      int new_bidx = wait_for_data(bundles);
      if (new_bidx < 0) break;
      bundle = &(bundles[new_bidx]);
#endif
    }

    // this_refid is an index into whatever array
    // refid is the actual refid for the specific read
    process_read_in(bundle->reads, id++, brec, strand, refid, 
      ref_name, brec->tag_int("NH"), brec->tag_int("HI"), hashread);

    prev_read_name = read_name;
  }
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
  std::string out_bam;
  std::string in_bam;
  std::string in_fasta;
  app.add_option("in.bam", in_bam, "input bam file")->required();
  // app.add_flag("--fr", FR_STRAND, "assume stranded library fw-firststrand");
  // app.add_flag("--rf", RF_STRAND, "assume stranded library fw-secondstrand");
  app.add_flag("--quiet", QUIET, "turn off verbose (log processing details)");
  app.add_flag("--long", LONG_READS, "alignments are from long reads?");
  app.add_option("-G", gff,
                 "reference annotation to use for guiding the assembly process "
                 "(GTF/GFF)")
      ->check(CLI::ExistingFile);
  app.add_option("-S", in_fasta,
                 "genome sequence for improving alignments "
                 "(FASTA)")
      ->check(CLI::ExistingFile);
  app.add_option("-o", out_bam,
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
  bam_file_in = in_bam.c_str();
  guide_gff = gff.c_str();
  bam_file_out = out_bam.c_str();

  if (!in_fasta.empty()) {
    gfasta = new GFastaDb(in_fasta.c_str());
    USE_FASTA = true;
  }

  // Add check for valid header

  // Create SAM header file
  const GStr header_path = "tx_header.sam";
  FILE *header_file = fopen(header_path.chars(), "w");
  if (header_file == NULL)
    GError("Error creating file: %s\n", header_path.chars());
  fprintf(header_file, "@HD\tVN:1.0\tSO:coordinate\n");

  if (!QUIET || BRAMBLE_DEBUG) {
    // print logo
    printf("\n-----------------------------------------------------\n");
    printf(" Bramble version %s \n", VERSION);
    printf("-----------------------------------------------------\n\n");
    printf("--> Loading reference annotation...\n");
  }
    
  // Open GFF/GTF file
  FILE *f = fopen(guide_gff.chars(), "r");
  if (f == NULL)
    GError("Error: could not open reference annotation file (%s)!\n",
           guide_gff.chars());

  // GffReader: transcripts only, sort by location
  GffReader *gffreader = new GffReader(f, true, true); // loading only recognizable transcript features
  gffreader->setRefAlphaSorted(); // alphabetical sorting of RefSeq IDs
  gffreader->showWarnings(BRAMBLE_DEBUG);

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
      if (BRAMBLE_DEBUG)
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

  fprintf(header_file, "@CO\tGenerated from GTF: %s\n", guide_gff.chars());
  fclose(header_file);

  if (!QUIET || BRAMBLE_DEBUG) {
    printf("Reference annotation loaded! %d unique transcripts were found\n\n", gffreader->gflst.Count());
  
    if (LONG_READS) {
      printf("--> Using LONG READS mode\n\n");
    }
    else {
      printf("--> Using SHORT READS mode\n");
      printf("Have long-read input data? Try running with the --long flag\n\n"); 
    }
    printf("--> Building g2t tree\n\n");
  }
  
  BamIO *io = new BamIO(bam_file_in, bam_file_out, header_path);
  io->start();

  auto g2t = build_g2t_tree(refguides, n_refguides, io);

  /*
  // Delete gffreader data
  for (int i = 0; i < gffreader->gflst.Count(); i++) {
    GffObj *guide = gffreader->gflst[i];
    delete guide;
  }
  */
  delete gffreader;

  #ifndef NOTHREADS

#define DEF_TSTACK_SIZE 8388608
  size_t def_stack_size = DEF_TSTACK_SIZE;

#ifdef _GTHREADS_POSIX_
  size_t tstack_size = GThread::defaultStackSize();
  if (tstack_size < DEF_TSTACK_SIZE)
    def_stack_size = DEF_TSTACK_SIZE;
  if (BRAMBLE_DEBUG) {
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

  if (!QUIET || BRAMBLE_DEBUG) {
    printf("--> Processing alignments\n\n");
  }

  process_reads(g2t, io, bundles, bundle, &bundle_queue);

   // Delete thread and bundle arrays
#ifndef NOTHREADS
  for (int t = 0; t < n_threads; t++) {
    threads[t].join();
  }
  if (BRAMBLE_DEBUG) {
    // LOG_INFO(logger, " All threads finished.\n");
  }
  delete[] threads;
  delete[] bundles;
#endif

  io->stop(); // close BAM reader & writer
  delete io;
  delete gfasta;

  printf("Final report:\n");
  printf("Total input reads:    %d\n", total_reads);
  printf("Unmapped reads:       %d\n", unmapped_reads);
  printf("Dropped reads:        %d\n", dropped_reads);
  printf("Unresolved reads:     %d\n\n", unresolved_reads);
  
  f_out = stdout;
  if (!QUIET || BRAMBLE_DEBUG)
    printf("-----------------------------------------------------\n");
}
