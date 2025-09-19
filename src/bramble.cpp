// Bramble v1.0.0

#include <numeric>

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
#include "bramble.h"
#include "reads.h"
#include "evaluate.h"

#define VERSION "1.0.0"

#define USAGE                                                                  \
  "Bramble v" VERSION " usage:\n\n\
bramble <in.bam ..> [-G <guide_gff>] [-o <out.bam>] [-p <cpus>] [-S <genome.fa>] \n\
 [--help] [--version] [--verbose] [--long] [--fr] [--rf] \n\
 \n\
A program to project spliced genomic alignments onto the transcriptome.\n\
Options:\n\
 --help     : print this usage message and exit\n\
 --version  : print just the version at stdout and exit\n\
 --verbose  : verbose (log bundle processing details)\n\
 --long     : BAM file contains long reads\n\
 --fr       : assume stranded library fw-firststrand\n\
 --rf       : assume stranded library fw-secondstrand\n\
 -G <file>  : reference annotation to use for guiding the BAM conversion (GTF/GFF)\n\
 -o <file>  : output path/file name for the projected alignments (default: stdout)\n\
 -p <int>   : number of threads (CPUs) to use (default: 1)\n\
 -S <file>  : genome sequence file (FASTA format)\n\
"

bool VERBOSE = false;     // Verbose, --verbose
bool LONG_READS = false;  // BAM file contains long reads, --long
bool FR_STRAND = true;   // Read 1 is on forward strand, --fr
bool RF_STRAND = false;   // Read 1 is on reverse strand, --fr
bool USE_FASTA = false;   // Use FASTA for reducing soft clips
bool SOFT_CLIPS = true;   // Add soft clips

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

// Check if there are more bundles
bool has_more_bundles() {
#ifndef NOTHREADS
    GLockGuard<GFastMutex> lock(reading_mutex);
#endif
  return !no_more_bundles;
}

// Wait for threads and set bundle status
void wait_for_bundles() {

#ifndef NOTHREADS
  reading_mutex.lock();
  no_more_bundles = true;
  reading_mutex.unlock();

  queue_mutex.lock();
  bundle_work &= ~(int)0x01; // clear bit 0;
  queue_mutex.unlock();

  bool are_threads_waiting = true;
  while (are_threads_waiting) {

    wait_mutex.lock();
    are_threads_waiting = (threads_waiting > 0);
    wait_mutex.unlock();

    if (are_threads_waiting) {
      have_bundles.notify_all();
      current_thread::sleep_for(1);

      wait_mutex.lock();
      are_threads_waiting = (threads_waiting > 0);
      wait_mutex.unlock();

      current_thread::sleep_for(1);
    }
  }
#else
  no_more_bundles = true;
#endif
}

// Process current bundle
void process_bundle(BundleData *bundle, BamIO *io) {

  if (VERBOSE) {
  #ifndef NOTHREADS
    GLockGuard<GFastMutex> lock(log_mutex);
  #endif
    GMessage(">bundle %s:%d-%d [%lu alignments (%d distinct), %d guides] begins processing...\n",
		        bundle->refseq.chars(), bundle->start, bundle->end, 
		        bundle->reads.Count(), bundle->guides.Count());
  
	}

  convert_reads(bundle, io);

  if (VERBOSE) {
    #ifndef NOTHREADS
    GLockGuard<GFastMutex> lock(log_mutex);
    #endif
    GMessage("^bundle %s:%d-%d done.\n",
             bundle->refseq.chars(), bundle->start, bundle->end);
  }

  bundle->Clear();
}

// Check that there aren't any threads waiting
bool no_threads_waiting() {
  wait_mutex.lock();
  int threads = threads_waiting;
  wait_mutex.unlock();
  return (threads < 1);
}

// Worker thread waits for incoming bundle, calls process_bundle
void worker_thread(GThreadData &td) {
  WorkerArgs *args = (WorkerArgs *)(td.udata);
  GPVec<BundleData> *bundle_queue = args->bundle_queue;
  BamIO *io = args->io;

  // Wait for a ready bundle in the queue
  queue_mutex.lock(); // enter wait-for-notification loop

  while (bundle_work) {
    wait_mutex.lock();
    threads_waiting++;
    queue_mutex.unlock();
    wait_mutex.unlock();
    have_threads.notify_one(); // in case main thread is waiting
    current_thread::yield();
    queue_mutex.lock();
    while (bundle_work && bundle_queue->Count() == 0) {
      // Unlocks queue_mutex and wait until notified
      // When notified, locks queue_mutex and resume
      have_bundles.wait(queue_mutex);
    }

    wait_mutex.lock();
    if (threads_waiting > 0)
      threads_waiting--;
    wait_mutex.unlock();

    BundleData *readyBundle = NULL;
    if ((bundle_work & 0x02) != 0) {
      readyBundle = bundle_queue->Pop();

      if (readyBundle != NULL) {
        if (bundle_queue->Count() == 0)
          bundle_work &= ~(int)0x02; // clear bit 1 (queue is empty)

        queue_mutex.unlock();
        process_bundle(readyBundle, io);

        data_mutex.lock();
        clear_data_pool.Push(readyBundle->idx);
        data_mutex.unlock();

        have_clear.notify_one(); // inform main thread
        current_thread::yield();

        queue_mutex.lock();
      }
    }
  }
  queue_mutex.unlock();
}

// Prepare the next available bundle slot for loading
int wait_for_data(BundleData *bundles) {
  int idx = -1;

  data_mutex.lock();
  while (clear_data_pool.Count() == 0) {
    have_clear.wait(data_mutex);
  }
  idx = clear_data_pool.Pop();
  if (idx >= 0)
    bundles[idx].status = BundleStatus::BUNDLE_STATUS_LOADING;
  data_mutex.unlock();

  return idx;
}

int main(int argc, char *argv[]) {

  //quill::Backend::start();
  // Get a logger
  // This segfaults on my system, not sure why
  // quill::Logger *logger = quill::Frontend::create_or_get_logger(
  //     "root",
  //     quill::Frontend::create_or_get_sink<quill::ConsoleSink>("sink_id_1"));

  CLI::App app{
      "A program to project spliced genomic alignments onto the transcriptome",
      "Bramble"};

  std::string gff;
  std::string bam_file;
  std::string input_bam;
  auto input_opt =
      app.add_option("in.bam", input_bam, "input bam file")->required();
  app.add_flag("--fr", FR_STRAND, "assume stranded library fw-firststrand");
  app.add_flag("--rf", RF_STRAND, "assume stranded library fw-secondstrand");
  app.add_flag("--verbose", VERBOSE, "verbose (log processing details)");
  app.add_flag("--long", LONG_READS, "verbose (log processing details)");
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
  // bam_file_in = new char[input_bam.size() + 1];
  // std::strcpy(bam_file_in, input_bam.c_str());
  bam_file_in = input_bam.c_str();
  guide_gff = gff.c_str();
  bam_file_out = bam_file.c_str();

  const char *ERR_BAM_SORT =
      "\nError: the input alignment file is not sorted!\n";

  // ^**^^^*^**^^^*^**^^^*^**^^^*
  // Read in reference annotation
  // ^**^^^*^**^^^*^**^^^*^**^^^*

  // Create SAM header file
  const GStr header_path = "tx_header.sam";
  FILE *header_file = fopen(header_path.chars(), "w");
  if (header_file == NULL)
    GError("Error creating file: %s\n", header_path.chars());
  fprintf(header_file, "@HD\tVN:1.0\tSO:coordinate\n");

  if (VERBOSE)
    GMessage(" Loading reference annotation (guides)..\n");

  // Open GFF/GTF file
  FILE *f = fopen(guide_gff.chars(), "r");
  if (f == NULL)
    GError("Error: could not open reference annotation file (%s)!\n",
           guide_gff.chars());

  // GffReader: transcripts only, sort by location
  GffReader *gffreader = new GffReader(f, true, true); // loading only recognizable transcript features
  gffreader->setRefAlphaSorted(); // alphabetical sorting of RefSeq IDs
  gffreader->showWarnings(VERBOSE);

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
      if (VERBOSE)
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

  if (VERBOSE) {
    GMessage("Transcriptome header created\n");
    GMessage(" %d reference transcripts loaded\n", gffreader->gflst.Count());
  }

  // ^**^^^*^**^^^*^**^^^*^**^^^*
  // Start BAM IO
  // ^**^^^*^**^^^*^**^^^*^**^^^*

  GMessage("long reads = %d\n", LONG_READS);

  BamIO *io = new BamIO(bam_file_in, bam_file_out, header_path);
  io->start();

  // *.** Bundle reads into overlapping groups

  GHash<int> hashread;          // read_name:pos:hit_index => readlist index
  GList<GffObj> *guides = NULL; // list of transcripts on a specific reference
  uint curr_bundle_start = 0;
  uint curr_bundle_end = 0;
  int first_possible_overlap = 0;
  int last_guide_idx = -1;
  int total_guides = 0;

  GStr last_ref_name;
  int last_ref_id = -1; // last seen ref_id

#ifndef NOTHREADS
#define DEF_TSTACK_SIZE 8388608
  size_t def_stack_size = DEF_TSTACK_SIZE;

#ifdef _GTHREADS_POSIX_
  size_t tstack_size = GThread::defaultStackSize();
  if (tstack_size < DEF_TSTACK_SIZE)
    def_stack_size = DEF_TSTACK_SIZE;
  if (VERBOSE) {
    if (tstack_size < def_stack_size) {
      GMessage("Default stack size for threads: %d (increased to %d)\n",
               tstack_size, def_stack_size);
    } else
      GMessage("Default stack size for threads: %d\n", tstack_size);
  }
#endif

  GThread *threads = new GThread[n_threads]; // threads for processing bundles
  GPVec<BundleData> *bundle_queue =
      new GPVec<BundleData>(false);          // queue for holding loaded bundles
  BundleData *bundles =
      new BundleData[n_threads + 1];         // redef with more bundles

  clear_data_pool.setCapacity(n_threads + 1);

  WorkerArgs *worker_args = new WorkerArgs(bundle_queue, io);

  // Start worker threads
  for (int b = 0; b < n_threads; b++) {
    threads[b].kickStart(worker_thread, (void *)worker_args, def_stack_size);
    bundles[b + 1].idx = b + 1;
    clear_data_pool.Push(b);
  }
  BundleData *bundle = &(bundles[n_threads]);

#else

  // Just put everything into the same bundle
  BundleData bundles[1];
  BundleData *bundle = &(bundles[0]);
#endif

  GSamRecord *brec = nullptr;
  bool more_alignments = true;
  int prev_pos = 0;

  while (more_alignments) {
    const char *ref_name = NULL;
    char splice_strand = 0;
    int read_start_pos = 0;
    int nh = 1;
    int hi = 0;
    int ref_id = last_ref_id; // current chromosome ID

    bool new_bundle = false;
    bool new_chromosome = false;


    if ((brec = io->next()) != NULL) {
      
      if (brec->isUnmapped())
        continue;

      ref_name = brec->refName();
      read_start_pos = brec->start;
      splice_strand = brec->spliceStrand(); // tagged strand gets priority

      // Set strand if stranded library
      if ((splice_strand == '.') && (FR_STRAND || RF_STRAND)) {

        if (brec->isPaired()) {
          if (brec->pairOrder() == 1) { // first read in pair
            if ((RF_STRAND && brec->revStrand()) ||
                (FR_STRAND && !brec->revStrand()))
              splice_strand = '+';
            else
              splice_strand = '-';
          } else { // second read in pair
            if ((RF_STRAND && brec->revStrand()) ||
                (FR_STRAND && !brec->revStrand()))
              splice_strand = '-';
            else
              splice_strand = '+';
          }

          // Read isn't paired
        } else {
          if ((RF_STRAND && brec->revStrand()) ||
              (FR_STRAND && !brec->revStrand()))
            splice_strand = '+';
          else
            splice_strand = '-';
        }
      }

      if (ref_name == NULL) {
        GError("Error: cannot retrieve target seq name from BAM record!\n");
      }

      // Are we at a new chromosome?
      new_chromosome = (last_ref_name.is_empty() || last_ref_name != ref_name);
      if (new_chromosome) {
        ref_id = guide_seq_names->gseqs.addName(ref_name);

        if (ref_id >= n_refguides) {
          if (VERBOSE) {
            GMessage("WARNING: no reference transcripts found for genomic "
                     "sequence \"%s\"!\n",
                     ref_name);
          }
        }

        prev_pos = 0;
      }
      if (read_start_pos < prev_pos)
        GError("%s\nread %s (start %d) found at position %d on %s when "
               "prev_pos=%d\n",
               ERR_BAM_SORT, brec->name(), brec->start, read_start_pos,
               ref_name, prev_pos);
      prev_pos = read_start_pos;

      nh = brec->tag_int("NH");
      nh = nh ? nh != 0 : 1;    // number of hits
      hi = brec->tag_int("HI"); // query hit index

      if ((!new_chromosome) && (curr_bundle_end > 0) &&
          (read_start_pos > (curr_bundle_end + (int)RUNOFF_DIST))) {
        new_bundle = true;
      }

      // No more alignments
    } else {
      more_alignments = false;
      new_bundle = true; // if last bundle, create fake start
    }

    // *.** New bundle / chromosome

    if (new_bundle || new_chromosome) {
      hashread.Clear();

      // Process reads in previous bundle
      if (bundle->reads.Count() > 0) {
        bundle->getReady(curr_bundle_start, curr_bundle_end);

        // For -S option: load genome FASTA
        if (USE_FASTA) {
          const char *chr_name = bundle->refseq.chars();
          GFaSeqGet *fasta_seq = gfasta->fetch(chr_name);

          // Try alternative naming convention
          if (fasta_seq == NULL) {
            GStr alt_chr_name;
            if (strncmp(chr_name, "chr", 3) == 0) {
              alt_chr_name = chr_name + 3; // skip "chr" prefix
            } else {
              alt_chr_name = "chr";
              alt_chr_name.append(chr_name); // add "chr" prefix
            }
            fasta_seq = gfasta->fetch(alt_chr_name.chars());
          }

          if (fasta_seq == NULL) {
            // LOG_ERROR(logger, "Error: could not retrieve sequence data for {}!\n",
            //           bundle->refseq.chars());
          }

          bundle->gseq = fasta_seq->copyRange(bundle->start, bundle->end, false, true);
        }

#ifndef NOTHREADS

        // Push this in the bundle queue where it'll be picked up by the threads
        int queue_count = 0;

        queue_mutex.lock();
        bundle_queue->Push(bundle);
        bundle_work |= 0x02; // set bit 1 to 1
        queue_count = bundle_queue->Count();
        queue_mutex.unlock();

        wait_mutex.lock();
        while (threads_waiting == 0) {
          have_threads.wait(wait_mutex);
        }
        wait_mutex.unlock();
        have_bundles.notify_one();

        current_thread::yield();

        queue_mutex.lock();
        while (bundle_queue->Count() == queue_count) {
          queue_mutex.unlock();
          have_bundles.notify_one();
          current_thread::yield();
          queue_mutex.lock();
        }
        queue_mutex.unlock();

#else // NOTHREADS = true

        // Just process single bundle
        process_bundle(bundle, io);
#endif
      }

      // Clear bundle (no more alignments)
      else {
        #ifndef NOTHREADS
          data_mutex.lock();
        #endif

        bundle->Clear();

        #ifndef NOTHREADS
          clear_data_pool.Push(bundle->idx);
          data_mutex.unlock();
        #endif
      }

      // If there is a new chromosome
      if (new_chromosome) {

        // Add guides from chromosome
        total_guides = 0;
        guides = nullptr;
        first_possible_overlap = 0;
        last_guide_idx = -1;

        // Check if BAM ref_id found in guide annotations
        if (refguides.Count() > ref_id) {

          // Check if refguides[ref_id] has guides
          if (refguides[ref_id].rnas.Count() > 0) {
            guides = &(refguides[ref_id].rnas);
            total_guides = guides->Count();
          } else {
            // LOG_WARNING(logger, "Warning: No guides for ref_id={} ({}), but ref present "
            //          "in guide index.\n",
            //          ref_id, ref_name);
          }
        } else {
          if (VERBOSE) {
            // LOG_WARNING(logger, "Warning: ref_id={} ({}) not found in guide annotations "
            //          "(refguides.Count() = {})\n",
            //          ref_id, ref_name, refguides.Count());
          }
        }

        last_ref_name = ref_name;
        last_ref_id = ref_id;
        curr_bundle_end = 0;
      }

      if (!more_alignments) {
        wait_for_bundles();
        break;
      }

      #ifndef NOTHREADS
        int new_bidx = wait_for_data(bundles);
        if (new_bidx < 0) {
          // LOG_ERROR(logger, "Error: wait_for_data() returned invalid bundle index({})!\n",
          //        new_bidx);
          break; // should never happen!
        }
        bundle = &(bundles[new_bidx]);
      #endif

      curr_bundle_start = read_start_pos;
      curr_bundle_end = brec->end;

      // *.** Add guides to bundle

      // Move to the first guide that could possibly overlap the read
      first_possible_overlap = last_guide_idx + 1;
      while (first_possible_overlap < total_guides &&
             static_cast<int>((*guides)[first_possible_overlap]->end) < read_start_pos) {
        first_possible_overlap++; // Skip non-overlapping guides
      }

      int guide_idx = first_possible_overlap;

      // Scan forward to find all guides overlapping the current read
      while (guide_idx < total_guides &&
             static_cast<int>((*guides)[guide_idx]->start) <= curr_bundle_end) {

        // Expand the current bundle range to include overlapping guides
        curr_bundle_start =
            std::min(curr_bundle_start, (*guides)[guide_idx]->start);
        curr_bundle_end = std::max(curr_bundle_end, (*guides)[guide_idx]->end);

        // On first overlapping guide, check for transitive overlaps *backward*
        if (guide_idx == first_possible_overlap && guide_idx > 0) {
          int back_idx = guide_idx;
          int first_transitive_overlap = guide_idx;

          // Look backward for any guide that also overlaps the updated region
          while (back_idx > last_guide_idx + 1) {
            --back_idx;
            if (curr_bundle_start <= (*guides)[back_idx]->end) {
              first_transitive_overlap = back_idx;
              curr_bundle_start =
                  std::min(curr_bundle_start, (*guides)[back_idx]->start);
            }
          }

          // Add all overlapping and transitively overlapping guides
          for (int i = first_transitive_overlap; i <= guide_idx; ++i) {
            bundle->keepGuide((*guides)[i]);
          }
        } else {
          // Directly overlapping guide
          bundle->keepGuide((*guides)[guide_idx]);
        }
        ++guide_idx;
      }

      // Update last processed guide index
      last_guide_idx = guide_idx - 1;

      bundle->refseq = last_ref_name;
      bundle->start = curr_bundle_start;
      bundle->end = curr_bundle_end;
    }

    // Current read extends the bundle
    if (curr_bundle_end < (int)brec->end) {
      curr_bundle_end = brec->end;

      // Add newly overlapping guides to the bundle
      // Keep expanding bundle to the right as long as new guides overlap with
      // the expanded region
      while (true) {
        bool expanded = false;
        while (last_guide_idx + 1 < total_guides &&
               (int)(*guides)[last_guide_idx + 1]->start <= curr_bundle_end) {
          ++last_guide_idx;
          auto *guide = (*guides)[last_guide_idx];
          if (guide->end >= curr_bundle_start) {
            bundle->keepGuide(guide);
            if (curr_bundle_end < guide->end) {
              curr_bundle_end = guide->end;
              expanded = true;
            }
          }
        }
        if (!expanded)
          break;
      }
    }

    process_read_in(curr_bundle_start, curr_bundle_end, *bundle, hashread, brec, splice_strand, nh, hi);
  }

  //^**^^^*^**^^^*^**^^^*^**^^^*
  // Finish and clean up
  //^**^^^*^**^^^*^**^^^*^**^^^*

  // Delete thread and bundle arrays
  #ifndef NOTHREADS
    for (int t = 0; t < n_threads; t++) {
      threads[t].join();
    }
    if (VERBOSE) {
      // LOG_INFO(logger, " All threads finished.\n");
    }
    delete[] threads;
    delete[] bundles;
    delete bundle_queue;
  #endif

  io->stop(); // close BAM reader & writer
  delete io;
  delete worker_args;
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
