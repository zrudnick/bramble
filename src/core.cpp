
#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include "types.h"
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
extern double similarity_threshold;

extern bool LONG_READS;
extern bool USE_FASTA;
extern bool SOFT_CLIPS;

namespace bramble {

  uint32_t get_mapq(uint32_t nh) {
    // MAPQ = int(-10log10(1-1/Nmap))
    if (!LONG_READS) {
      if (nh == 1) return 255;
      else if (nh == 2) return 3;
      else if (nh == 3 || nh == 4) return 1;
      else return 0;
    // LONG READS: mapq > 1 if there are secondary alignments
    } else {
      if (nh > 1) return 0;
      else return 3;
    }
  }

  ReadInfo *process_read_out(CReadAln &read, read_id_t id, 
                            std::shared_ptr<g2tTree> g2t, 
                            std::shared_ptr<ReadEvaluator> evaluator,
                            uint8_t *seq, int seq_len) {

    // Evaluate read group to find matches
    std::vector<ExonChainMatch> matches = 
      evaluator->evaluate(read, g2t, seq, seq_len);

    if (matches.empty()) return nullptr;

    // Write out evaluation result
    auto this_read = new ReadInfo;

    this_read->matches = std::move(matches);
    this_read->valid_read = true;
    this_read->is_paired = (read.brec->flags() & BAM_FPAIRED);
    this_read->read = std::make_shared<ReadOut>();
    this_read->read->index = id; // used for seen in bam
    this_read->read->brec = read.brec;
    this_read->read->nh = matches.size();
    return this_read;
  }

  void write_to_bam(BamIO *io, std::vector<BamInfo *>& bam_info) {
    unordered_set<read_id_t> seen;
    CigarMem cigar_mem;             // reuse space in memory
    std::vector<bam1_t*> to_write;  // collect finished records

    // Prepare records for paired reads or single read
    for (const auto& pair : bam_info) {
      if (!pair || !pair->valid_pair) continue;

      // Prepare records for one output alignment
      auto prepare_read = [&](std::shared_ptr<ReadOut> read, 
                              Cigar &cigar,
                              char strand,
                              bool is_first) {
        if (!read || !read->brec) return;

        bam1_t* old_b = read->brec->get_b();
        bam1_t* b; // new copy

        bool new_read = seen.insert(read->index).second;
        if (new_read && !LONG_READS) {
          // Set tags
          set_nh_tag(old_b, read->nh);
          set_xs_tag(old_b, strand);
        } else if (new_read && LONG_READS) { 
          // Set NH and XS tags
          set_nh_tag(old_b, read->nh);
          set_ts_tag(old_b, strand);
        }

        b = bam_dup1(old_b);

        // Print BEFORE update
        uint32_t* new_cigar = bam_get_cigar(b);
        uint32_t n_cigar = b->core.n_cigar;
        int32_t nm = 0;
        update_cigar(b, new_cigar, n_cigar, cigar_mem, cigar, nm);
        // set_nm_tag(b, nm); should this be updated?

        // Update mapping quality score
        // must have values between 0–255
        b->core.qual = (uint8_t)read->mapq;
           
        // Set coordinates and align-specific information
        if (is_first) {
          b->core.tid = (int32_t)pair->r_tid;
          if (pair->r_align.primary_alignment) b->core.flag &= ~BAM_FSECONDARY; // mark primary
          else b->core.flag |= BAM_FSECONDARY; // default secondary
          // If read is being reversed - reverse-complement the sequence and quality strings
          if (strand == '-') {
            // update to allow * for sequence
            int rc_res = reverse_complement_bam(b);
            if (rc_res != 0) {
              GError("Error: reverse_complement_bam failed with code %d\n", rc_res);
            }
          }

          if (strand == '+') { // pos depends only on transcript strand
                               // not on read orientation or reverse status
            b->core.pos = (int32_t)pair->r_align.fwpos;
          } else {
            b->core.pos = (int32_t)pair->r_align.rcpos;
          }
          
          if (LONG_READS) set_as_tag(b, pair->r_align, b->core.qual);
          set_hi_tag(b, pair->r_align.hit_index);

        } else {  // mate
          b->core.tid = (int32_t)pair->m_tid;
          if (pair->m_align.primary_alignment) b->core.flag &= ~BAM_FSECONDARY; // mark primary
          else b->core.flag |= BAM_FSECONDARY; // default secondary
          // If read is being reversed - reverse-complement the sequence and quality strings
          if (strand == '-') {
            // update to allow * for sequence
            int rc_res = reverse_complement_bam(b);
            if (rc_res != 0) {
              GError("Error: reverse_complement_bam failed with code %d\n", rc_res);
            }
          }

          if (strand == '+') { // pos depends only on transcript strand
                               // not on read orientation or reverse status
            b->core.pos = (int32_t)pair->m_align.fwpos;
          } else {
            b->core.pos = (int32_t)pair->m_align.rcpos;
          }
          
          if (LONG_READS) set_as_tag(b, pair->m_align, b->core.qual);
          set_hi_tag(b, pair->m_align.hit_index);
        }

        set_mate_info(b, pair, is_first); // must be called for both paired and non-paired instances
        to_write.push_back(b);

      };

      prepare_read(pair->read1, pair->r_align.cigar, 
        pair->r_align.strand, true);
      if (pair->is_paired)
        prepare_read(pair->read2, pair->m_align.cigar, 
          pair->m_align.strand, false);
    }

    // Write records together to enable reasonable BAM compression
#ifndef NOTHREADS
    bam_io_mutex.lock();
#endif

    for (auto* b_ : to_write) {
      io->write(b_);
      bam_destroy1(b_);
    }

#ifndef NOTHREADS
    bam_io_mutex.unlock();
#endif
  }

  void convert_reads(std::vector<CReadAln> &reads,
                    std::shared_ptr<g2tTree> g2t, 
                    std::shared_ptr<ReadEvaluator> evaluator, 
                    BamIO *io) {

    // Chunk buffer for BAM output
    std::vector<BamInfo*> bam_info;
    bam_info.reserve(CHUNK_SIZE*1.2);

    // Buffer for grouping by read name before filtering
    std::unordered_map<const char *, std::vector<BamInfo*>> pairs_by_name;
    uint32_t n_pairs = 0;

    auto flush = [&]() {
      // Filter by strand before writing
      std::vector<BamInfo*> filtered_bam_info;
      
      for (auto& [read_name, pairs] : pairs_by_name) {
      
        // Recalculate NH and MAPQ based on number of kept alignments
        uint32_t new_nh = pairs.size();
        uint32_t new_mapq = get_mapq(new_nh);
        
        for (auto* pair : pairs) {
          pair->read1->nh = new_nh;
          pair->read1->mapq = new_mapq;
          filtered_bam_info.push_back(pair);
        }
      }
      
      write_to_bam(io, filtered_bam_info);
      for (auto* info : filtered_bam_info) delete info;
      pairs_by_name.clear();
      n_pairs = 0;
    };

    auto emit_pair = [&](BamInfo* pair, bool is_last) {
      if (pair && pair->read1) {
        const char *read_name = pair->read1->brec->name();
        pairs_by_name[read_name].push_back(pair);
        n_pairs++;
        
        if (is_last && n_pairs >= CHUNK_SIZE) {
          flush();
        }
      }
    };

    // Track seen reads
    unordered_set<read_id_t> seen;
    seen.reserve(reads.size() * 1.2);

    for (read_id_t id = 0; id < reads.size(); ) {

      // Group reads with identical read names
      read_id_t start = id;
      const std::string& name = reads[id].brec->name();
      id++; // prevent doing it twice

      while (id < reads.size() && reads[id].brec->name() == name) {
        id++;
      }
      read_id_t end = id; // [start, end)

      CReadAln &first = reads[start];

      int seq_len;
      uint8_t *seq;
      if (LONG_READS) {
        bam1_t* b = first.brec->get_b();
        seq_len = b->core.l_qseq;
        seq = bam_get_seq(b);
      } else {
        seq_len = 0;
        seq = nullptr;
      }
      
      // Process read groups
      for (read_id_t i = start; i < end; i++) {
        if (seen.count(i)) continue;

        ReadInfo* this_read = process_read_out(
          reads[i], i, g2t, evaluator, seq, seq_len);

        int n_mates = reads[i].pair_idx.Count();
        if (n_mates == 0) {
          process_mate_pair(this_read, nullptr, emit_pair);
          delete this_read;
          seen.insert(i);
          continue;
        }

        for (uint32_t m = 0; m < n_mates; m++) {
          int mate_id = reads[i].pair_idx[m];
          if (mate_id < 0 || mate_id >= reads.size()) continue;
          if (seen.count(mate_id)) continue;

          ReadInfo* mate_read = process_read_out(
            reads[mate_id], mate_id, g2t, evaluator, seq, seq_len);

          process_mate_pair(this_read, mate_read, emit_pair);
          delete mate_read;
          seen.insert(mate_id);
        }

        delete this_read;
        seen.insert(i);
      }
    }

    if (!pairs_by_name.empty()) flush();
  }

}