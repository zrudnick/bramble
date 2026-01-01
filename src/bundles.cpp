
#include <numeric>
#include <algorithm>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <iomanip>
#include <string>
#include <random>
#include <chrono>

#ifndef NOTHREADS
#include "GThreads.h"
#endif
#include "htslib/sam.h"

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

extern bool VERBOSE;     // Verbose, --verbose
extern bool DEBUG_;
extern bool LONG_READS;  // BAM file contains long reads, --long
extern bool FR_STRAND;   // Read 1 is on forward strand, --fr
extern bool RF_STRAND;   // Read 1 is on reverse strand, --fr
extern bool USE_FASTA;   // Use FASTA for reducing soft clips
extern bool SOFT_CLIPS;   // Add soft clips

const int MAX_READS_PER_BUNDLE = 50000;

extern uint8_t n_threads;    // Threads, -p
extern GFastaDb *gfasta;  // FASTA file, -S

#ifndef NOTHREADS
// Threading: single producer, multiple consumers
// Main thread/program is always loading the producer

extern GMutex data_mutex; // Manage availability of data records ready to be loaded by
                   // main thread
extern GVec<int> clear_data_pool; // Indices of data bundles cleared for loading by
                           // main thread (clear data pool)
extern GConditionVar
    have_bundles; // Will notify a thread that a bundle was loaded in the ready
                  // queue (or that no more bundles are coming)
extern char bundle_work; // Bit 0 set if bundles are still being prepared (BAM file not exhausted
       // yet) Bit 1 set if there are Bundles ready in the queue

extern GMutex wait_mutex;       // Controls threads_waiting (idle threads counter)
extern uint8_t threads_waiting; // Idle worker threads
extern GConditionVar
    have_threads; // Will notify the bundle loader when a thread
                  // Is available to process the currently loaded bundle

extern GConditionVar have_clear; // Will notify when bundle buf space available
extern GMutex queue_mutex;       // Controls bundle_queue and bundles access
extern GFastMutex log_mutex;     // Only when verbose - to avoid mangling the log output
extern GFastMutex reading_mutex;
extern GFastMutex bam_io_mutex;  // Protects BAM io
#endif

extern bool no_more_bundles;

// uint32_t dropped_reads;
// uint32_t unresolved_reads;

namespace bramble {

  void process_exons(GSamRecord *brec, CReadAln *read) {
    auto exons = brec->exons;
    for (int i = 0; i < exons.Count(); i++) {
      read->len += exons[i].len();
      read->segs.Add(exons[i]);
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
                      int32_t nh, int32_t hi, GHash<int> &hashread) {

    // Create new read alignment
    CReadAln *read = new CReadAln(strand, refid, nh, brec->start, brec->end);
    read->brec = brec;
    reads.emplace_back(read);
    
    process_exons(brec, read);
    float read_count = calculate_read_count(brec);
    process_pairs(reads, id, brec, hi, read_count, hashread);
  }

  char get_splice_strand(GSamRecord* brec){
    char splice_strand = brec->spliceStrand(); // tagged strand gets priority

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

    return splice_strand;
  }

  // Process reads in previous bundle
  void push_bundle(BundleData *bundle, GPVec<BundleData> *bundle_queue,
                  int curr_bundle_start, int curr_bundle_end) {

    if (bundle->reads.size() > 0) {
      bundle->getReady(curr_bundle_start, curr_bundle_end);

      // For -S option: load genome FASTA
      if (USE_FASTA) {
        const char *chr_name = bundle->refseq.chars();
        GFaSeqGet *fasta_seq = gfasta->fetch(chr_name);

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
  }

  void build_bundles(BamIO* io, BundleData *bundles, BundleData* bundle, 
                    GPVec<BundleData> *bundle_queue, GVec<GRefData> refguides, 
                    int n_refguides, std::shared_ptr<GffNames> guide_seq_names) {

    GHash<int> hashread;
    GList<GffObj> *guides = NULL;
    int total_guides = 0;

    GStr last_ref_name;
    int last_ref_id = -1;

    GSamRecord *brec = nullptr;
    int rec_count = 0;
    
    bool more_alignments = true;
    int prev_pos = 0;

    // Track the bundle boundaries
    uint bundle_min_pos = 0;
    uint bundle_max_pos = 0;
    
    // Track which guides we've already added to avoid duplicates
    int first_guide_in_range = 0;

    // uint32_t all_reads = 0;
    // uint32_t unmapped_reads = 0;
    // dropped_reads = 0;
    // unresolved_reads = 0;

    while (more_alignments) {
      const char *ref_name = NULL;
      char splice_strand = '.';
      int read_start_pos = 0;
      int nh = 1;
      int hi = 0;
      int ref_id = last_ref_id;

      bool new_bundle = false;
      bool new_chromosome = false;

      // Get next alignment
      if ((brec = io->next()) != NULL) {
        rec_count += 1;
        //all_reads++;
        
        if (brec->isUnmapped()) {
          //unmapped_reads++;
          continue;
        }

        ref_name = brec->refName();
        read_start_pos = brec->start;
        splice_strand = get_splice_strand(brec);

        if (ref_name == NULL) {
          GError("Error: cannot retrieve target seq name from BAM record!\n");
        }

        // Check for chromosome change
        new_chromosome = (last_ref_name.is_empty() || last_ref_name != ref_name);
        if (new_chromosome) {
          ref_id = guide_seq_names->gseqs.addName(ref_name);
          if (ref_id >= n_refguides && DEBUG_) {
            GMessage("WARNING: no reference transcripts found for genomic sequence \"%s\"!\n", ref_name);
          }
          prev_pos = 0;
        }

        // Validate sorted order
        if (read_start_pos < prev_pos) {
          GError("Error: the input alignment file is not sorted!\n"
                 "read %s (start %d) found at position %d on %s when prev_pos=%d\n",
                 brec->name(), brec->start, read_start_pos, ref_name, prev_pos);
        }
        prev_pos = read_start_pos;

        nh = brec->tag_int("NH");
        nh = nh ? nh : 1;
        hi = brec->tag_int("HI");

        // Check if we need a new bundle
        if (!new_chromosome && bundle_max_pos > 0 && 
            read_start_pos > bundle_max_pos + RUNOFF_DIST) {
          new_bundle = true;
        }

        // Check if we hit max reads per bundle
        if (rec_count >= MAX_READS_PER_BUNDLE) {
          new_bundle = true;
        }

      } else {
        // No more alignments
        more_alignments = false;
        new_bundle = true;
      }

      // Process bundle boundary
      if (new_bundle || new_chromosome) {
        hashread.Clear();
        rec_count = 0;

        // Push completed bundle
        push_bundle(bundle, bundle_queue, bundle_min_pos, bundle_max_pos);

        if (new_chromosome) {
          total_guides = 0;
          guides = nullptr;
          first_guide_in_range = 0;

          // Load guides for new chromosome
          if (refguides.Count() > ref_id && refguides[ref_id].rnas.Count() > 0) {
            guides = &(refguides[ref_id].rnas);
            total_guides = guides->Count();
          }

          last_ref_name = ref_name;
          last_ref_id = ref_id;
        }

        if (!more_alignments) {
          //wait_for_bundles();
          break;
        }

        #ifndef NOTHREADS
          //int new_bidx = wait_for_data(bundles);
          //if (new_bidx < 0) break;
          //bundle = &(bundles[new_bidx]);
        #endif

        // Initialize new bundle with first read
        bundle_min_pos = read_start_pos;
        bundle_max_pos = brec->end + 1; // exclusive

        // Find guides that overlap this bundle
        // Since reads are sorted and we broke the bundle, we can skip guides that ended before this position
        if (guides != nullptr) {
          // Move first_guide_in_range forward to skip guides that definitely can't overlap
          while (first_guide_in_range < total_guides && 
                 (*guides)[first_guide_in_range]->end <= bundle_min_pos) {
            first_guide_in_range++;
          }
          
          // Now scan from first_guide_in_range and add overlapping guides
          // Keep expanding until no new guides are found
          bool expanded = true;
          while (expanded) {
            expanded = false;
            int scan_idx = first_guide_in_range;
            
            while (scan_idx < total_guides) {
              uint guide_start = (*guides)[scan_idx]->start;
              uint guide_end = (*guides)[scan_idx]->end + 1;
              
              // If guide starts after bundle ends, we're done (guides are sorted)
              if (guide_start >= bundle_max_pos) {
                break;
              }
              
              // Check if guide overlaps current bundle range
              if (guide_end > bundle_min_pos) {
                bundle->keepGuide((*guides)[scan_idx]);
                
                // Expand bundle to include this guide
                if (guide_start < bundle_min_pos) {
                  bundle_min_pos = guide_start;
                  expanded = true;
                }
                if (guide_end > bundle_max_pos) {
                  bundle_max_pos = guide_end;
                  expanded = true;
                }
              }
              
              scan_idx++;
            }
          }
        }

        bundle->refseq = last_ref_name;
        bundle->start = bundle_min_pos;
        bundle->end = bundle_max_pos;
      }

      // Extend bundle with current read
      if (brec->end + 1 > bundle_max_pos) {
        uint old_max = bundle_max_pos;
        bundle_max_pos = brec->end + 1;

        // Add any newly overlapping guides
        // Only need to check guides that start between old_max and new bundle_max_pos
        if (guides != nullptr) {
          bool expanded = true;
          while (expanded) {
            expanded = false;
            int scan_idx = first_guide_in_range;
            
            while (scan_idx < total_guides) {
              uint guide_start = (*guides)[scan_idx]->start;
              uint guide_end = (*guides)[scan_idx]->end + 1;
              
              // If guide starts after bundle ends, we're done
              if (guide_start >= bundle_max_pos) {
                break;
              }
              
              // Check if guide overlaps current bundle range
              if (guide_end > bundle_min_pos) {
                bundle->keepGuide((*guides)[scan_idx]);
                
                if (guide_end > bundle_max_pos) {
                  bundle_max_pos = guide_end;
                  expanded = true;
                }
              }
              
              scan_idx++;
            }
          }
        }
      }

      // Process the read
      // process_read_in(bundle_min_pos, bundle_max_pos, *bundle, 
      //   hashread, brec, splice_strand, nh, hi);
    }

    // GMessage("#total input reads is:    %d\n", all_reads);
    // GMessage("#unmapped reads is: %d\n", unmapped_reads);
    // GMessage("#dropped reads is: %d\n", dropped_reads);
    // GMessage("#unresolved reads is: %d\n\n", unresolved_reads);
  }

}