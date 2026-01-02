
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
const uint32_t CHUNK_SIZE = 50000;   // size of each write to BAM

extern bool LONG_READS;
extern bool USE_FASTA;
extern bool SOFT_CLIPS;

namespace bramble {

  // MAPQ = int(-10log10(1-1/Nmap))
  uint32_t get_mapq(uint32_t nh) {
    if (nh == 1) return 255;
    else if (nh == 2) return 3;
    else if (nh == 3 || nh == 4) return 1;
    else return 0;
  }

  ReadInfo *process_read_out(CReadAln *read, read_id_t id, 
                            g2tTree *g2t, ReadEvaluator *evaluator) {

    // Evaluate read group to find matches
    std::unordered_map<tid_t, ExonChainMatch> matches = 
      evaluator->evaluate(read, id, g2t);

    if (matches.empty()) return nullptr;

    // Write out evaluation result
    auto this_read = new ReadInfo;

    this_read->matches = matches;
    this_read->valid_read = true;
    this_read->is_paired = (read->brec->flags() & BAM_FPAIRED);
    this_read->read = std::make_shared<ReadOut>();
    this_read->read->name = read->brec->name();
    this_read->read->index = id;
    this_read->read->brec = read->brec;
    this_read->read->read_size = read->len;
    this_read->read->nh = matches.size();
    this_read->read->mapq = get_mapq(this_read->read->nh);
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
        // if (new_read && !LONG_READS) {
        //   uint32_t* cigar = bam_get_cigar(old_b);
        //   uint32_t n_cigar = old_b->core.n_cigar;

        //   int32_t nm = 0;
        //   update_cigar(old_b, cigar, n_cigar, cigar_mem, *ideal_cigar, nm);

        //   // Set tags
        //   set_nh_tag(old_b, read->nh);
        //   set_xs_tag(old_b, strand);
        //   set_nm_tag(old_b, nm);
        //   remove_extra_tags(old_b);

        // we have different cigars for every match
        // so just modify new b every time
        // } else if (new_read && LONG_READS) { 
          // Set NH and XS tags
          set_nh_tag(old_b, read->nh);
          set_xs_tag(old_b, strand);
          remove_extra_tags(old_b);
        // }

        // if (!LONG_READS) {
        //   b = bam_dup1(old_b);
        // } else {

          b = bam_dup1(old_b);

          // Print BEFORE update
          uint32_t* cigar = bam_get_cigar(b);
          uint32_t n_cigar = b->core.n_cigar;
          int32_t nm = 0;
          update_cigar(b, cigar, n_cigar, cigar_mem, *ideal_cigar, nm);
          set_nm_tag(b, nm);
        // }
           
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

        // Update mapping quality score
        // must have values between 0–255
        uint8_t prev_qual = b->core.qual;
        if (prev_qual <= 3) {
          if (read->mapq > 3) b->core.qual = prev_qual;
          else b->core.qual = (uint8_t)read->mapq;
        } else {
          b->core.qual = (uint8_t)read->mapq;
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

  void convert_reads(std::vector<CReadAln *> &reads,
                    g2tTree* g2t, ReadEvaluator* evaluator, 
                    BamIO *io) {

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

    for (read_id_t id = 0; id < reads.size(); id++) {
      if (seen.count(id)) continue; // already processed as a mate

      ReadInfo* this_read = process_read_out(reads[id], 
        id, g2t, evaluator);
      int n_mates = reads[id]->pair_idx.Count();

      // Unpaired reads
      if (n_mates == 0) {
        process_mate_pair(this_read, nullptr, emit_pair);
        delete this_read;
        seen.insert(id);
        continue;
      }

      // Paired reads: process this_read with all its mates
      for (uint32_t m = 0; m < n_mates; m++) {
        int mate_id = reads[id]->pair_idx[m];
        if (mate_id < 0 || mate_id >= reads.size()) continue;

        if (seen.count(mate_id)) continue;

        ReadInfo* mate_read = process_read_out(reads[mate_id], 
          mate_id, g2t, evaluator);
        process_mate_pair(this_read, mate_read, emit_pair);
        delete mate_read;
        seen.insert(mate_id);
      }

      delete this_read;
      seen.insert(id);
    }

    if (!bam_info.empty()) flush();
  }

}