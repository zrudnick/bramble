
#include <algorithm>
#include <memory>
#include <set>
#include <string>
#include <vector>
#include <unordered_set>

#include "types.h"
#include "bundles.h"
#include "bramble.h"
#include "g2t.h"
#include "evaluate.h"
#include "bam.h"
#include "mates.h"

#ifndef NOTHREADS
#include "GThreads.h"
#endif

extern GFastMutex bam_io_mutex;     // protects BAM io
const uint32_t CHUNK_SIZE = 5000;   // size of each write to BAM

extern bool LONG_READS;
extern bool USE_FASTA;
extern bool SOFT_CLIPS;

namespace bramble {
  /**
   * Prints g2t tree with nodes sorted by genome coordinate
   * For debugging only
   *
   * @param g2t g2t tree for bundle
   */
  void print_tree(g2tTree* g2t) {
    auto intervals = g2t->fw_tree->getOrderedIntervals();
    printf("FW Tree:\n");
    for (auto &interval : intervals) {
      printf("(%u, %u)", interval->start, interval->end);
      printf(" size = %u, %u", interval->end - interval->start);
      for (auto &tid : interval->tids) {
        printf(" %s ", g2t->getTidName(tid).c_str());
        printf(" %d ", interval->tid_cum_len[tid]);
      }
      printf("-> \n");
    }
    printf("END\n");

    intervals = g2t->rc_tree->getOrderedIntervals();
    printf("RC Tree:\n");
    for (auto &interval : intervals) {
      printf("(%u, %u)", interval->start, interval->end);
      for (auto &tid : interval->tids) {
        printf(" %s ", g2t->getTidName(tid).c_str());
        printf(" %d ", interval->tid_cum_len[tid]);
      }
      printf("-> \n");
    }
    printf("END\n");
  }

  /**
   * Builds g2t tree from guide transcripts
   *
   * @param bundle current bundle of reads
   * @return complete g2t tree for bundle
   */
  std::unique_ptr<g2tTree> make_g2t_tree(BundleData *bundle, BamIO *io) {
    GPVec<GffObj> guides = bundle->guides;

    if (guides.Count() == 0)
      return nullptr;

    auto g2t = std::make_unique<g2tTree>();

    // Insert all guide exons into the interval tree
    for (int i = 0; i < guides.Count(); i++) {
      GffObj *guide = guides[i];
      char strand = guide->strand;
      const char *tid_string = guide->getID(); // we don't strictly need this

#ifndef NOTHREADS
      bam_io_mutex.lock();
#endif
      tid_t tid = g2t->insertTidString(tid_string, io);
#ifndef NOTHREADS
      bam_io_mutex.unlock();
#endif

      for (int j = 0; j < guide->exons.Count(); j++) {
        uint g_start, g_end;
        if (strand == '+') {
          g_start = guide->exons[j]->start - bundle->start;
          g_end = (guide->exons[j]->end + 1) - bundle->start;
        } else {
          g_start = bundle->end - (guide->exons[j]->end + 1);
          g_end = bundle->end - guide->exons[j]->start;
        }

        g2t->addGuideExon(g_start, g_end, tid, strand);
      }
    }

    g2t->buildAllTidChains();
    g2t->precomputeAllCumulativeLengths();

// #ifndef NOTHREADS
//     bam_io_mutex.lock();
// #endif
//     print_tree(g2t.get());
// #ifndef NOTHREADS
//     bam_io_mutex.unlock();
// #endif
    return g2t;
  }

  ReadInfo *process_read_out(BundleData *&bundle, uint read_idx, g2tTree *g2t,
                            ReadEvaluator *evaluator) {

    // Evaluate read group to find matches
    ReadEvaluationResult res = evaluator->evaluate(bundle, read_idx, g2t);

    if (!res.valid) return nullptr;

    // Write out evaluation result
    auto read = bundle->reads[read_idx];
    auto this_read = new ReadInfo;

    this_read->matches = res.matches;
    this_read->valid_read = true;
    this_read->is_paired = (read->brec->flags() & BAM_FPAIRED);
    this_read->read = std::make_shared<ReadOut>();
    this_read->read->name = read->brec->name();
    this_read->read->index = read_idx;
    this_read->read->brec = read->brec;
    this_read->read->read_size = read->len;
    this_read->read->nh = res.matches.size();

    return this_read;
  }

  void write_to_bam(BamIO *io, std::vector<BamInfo *>& bam_info) {
    std::unordered_set<read_id_t> seen;
    CigarMem cigar_mem; // reuse space in memory
    std::vector<bam1_t*> to_write;  // collect finished records

    // Prepare records for paired reads or single read
    for (const auto& this_pair : bam_info) {
      if (!this_pair || !this_pair->valid_pair) continue;

      // Prepare records for one output alignment
      auto prepare_read = [&](std::shared_ptr<ReadOut> read, 
                              std::shared_ptr<Cigar> &ideal_cigar,
                              char strand,
                              bool is_first) {
        if (!read || !read->brec) return;

        bam1_t* old_b = read->brec->get_b();
        bam1_t* b; // new copy

        bool new_read = seen.insert(read->index).second;
        if (new_read && !LONG_READS) {
          uint32_t* cigar = bam_get_cigar(old_b);
          uint32_t n_cigar = old_b->core.n_cigar;

          int32_t nm = 0;
          update_cigar(old_b, cigar, n_cigar, cigar_mem, *ideal_cigar, nm);

          // Set tags
          set_nh_tag(old_b, read->nh);
          set_xs_tag(old_b, strand);
          set_nm_tag(old_b, nm);
          remove_extra_tags(old_b);

        // we have different cigars for every match
        // so just modify new b every time
        } else if (new_read && LONG_READS) { 
          // Set NH and XS tags
          set_nh_tag(old_b, read->nh);
          set_xs_tag(old_b, strand);
          remove_extra_tags(old_b);
        }

        if (!LONG_READS) {
          b = bam_dup1(old_b);
        } else {

          b = bam_dup1(old_b);

          // Print BEFORE update
          uint32_t* cigar = bam_get_cigar(b);
          uint32_t n_cigar = b->core.n_cigar;
          int32_t nm = 0;
          update_cigar(b, cigar, n_cigar, cigar_mem, *ideal_cigar, nm);
          set_nm_tag(b, nm);
        }
           
        // Set coordinates and align-secific information
        if (is_first) {
          b->core.tid = (int32_t)this_pair->r_tid;
          b->core.pos = (int32_t)this_pair->r_align.pos;
          if (this_pair->r_align.primary_alignment) b->core.flag &= ~BAM_FSECONDARY; // mark primary
          else b->core.flag |= BAM_FSECONDARY; // default secondary
          
          set_as_tag(b, this_pair->r_align.similarity_score);
          set_hi_tag(b, this_pair->r_align.hit_index);
        } else {
          b->core.tid = (int32_t)this_pair->m_tid;
          b->core.pos = (int32_t)this_pair->m_align.pos;
          if (this_pair->m_align.primary_alignment) b->core.flag &= ~BAM_FSECONDARY; // mark primary
          else b->core.flag |= BAM_FSECONDARY; // default secondary
          
          set_as_tag(b, this_pair->m_align.similarity_score);
          set_hi_tag(b, this_pair->m_align.hit_index);
        }

        // If read is being reversed - reverse-complement the sequence and quality strings
        if (b->core.flag & BAM_FREVERSE) {
          // update to allow * for sequence
          int rc_res = reverse_complement_bam(b);
          if (rc_res != 0) {
            GError("Error: reverse_complement_bam failed with code %d\n", rc_res);
          }
        }

        set_mate_info(b, this_pair, is_first); // must be called for both paired and non-paired instances

        to_write.push_back(b);

      };

      prepare_read(this_pair->read1, this_pair->r_align.cigar, 
        this_pair->r_align.strand, true);
      if (this_pair->is_paired)
        prepare_read(this_pair->read2, this_pair->m_align.cigar, 
          this_pair->m_align.strand, false);
    }

    // Write records together to enable reasonable BAM compression
#ifndef NOTHREADS
    bam_io_mutex.lock();
#endif

    for (auto* b_ : to_write) {
      GSamRecord *rec = new GSamRecord(b_);
      io->write(rec);
      delete rec;
    }

#ifndef NOTHREADS
    bam_io_mutex.unlock();
#endif
  }

  void convert_reads(BundleData *bundle, BamIO *io) {
    // Build g2t tree
    std::unique_ptr<g2tTree> g2t = make_g2t_tree(bundle, io);
    if (g2t == nullptr) return;

    // Pick evaluator
    std::unique_ptr<ReadEvaluator> evaluator;
    if (LONG_READS) {
      evaluator = std::make_unique<LongReadEvaluator>();
    } else {
      evaluator = std::make_unique<ShortReadEvaluator>();
    }

    std::vector<CReadAln *> &reads = bundle->reads;

    // Chunk buffer for BAM output
    std::vector<BamInfo*> bam_info;
    bam_info.reserve(CHUNK_SIZE*1.2);

    auto flush = [&]() {
      write_to_bam(io, bam_info);
      for (auto* info : bam_info) delete info;
      bam_info.clear();
    };

    auto emit_pair = [&](BamInfo* pair, bool is_last) {
      bam_info.push_back(pair);
      if (is_last && bam_info.size() >= CHUNK_SIZE) {
        flush();
      }
    };

    // Track seen reads
    std::unordered_set<read_id_t> seen;
    seen.reserve(reads.size() * 1.2);

    for (read_id_t n = 0; n < reads.size(); n++) {
      if (seen.count(n)) continue; // already processed as a mate

      ReadInfo* this_read = process_read_out(bundle, 
        n, g2t.get(), evaluator.get());
      int n_mates = reads[n]->pair_idx.Count();

      // Unpaired reads
      if (n_mates == 0) {
        process_mate_pair(this_read, nullptr, emit_pair);
        delete this_read;
        seen.insert(n);
        continue;
      }

      // Paired reads: process this_read with all its mates
      for (uint32_t m = 0; m < n_mates; m++) {
        int mate_idx = reads[n]->pair_idx[m];
        if (mate_idx < 0 || mate_idx >= reads.size()) continue;

        if (seen.count(mate_idx)) continue;

        ReadInfo* mate_read = process_read_out(bundle, 
          mate_idx, g2t.get(), evaluator.get());
        process_mate_pair(this_read, mate_read, emit_pair);
        delete mate_read;
        seen.insert(mate_idx);
      }

      delete this_read;
      seen.insert(n);
    }

    if (!bam_info.empty()) flush();
  }

}