
#include <algorithm>
#include <memory>
#include <string>
#include <vector>
#include <cmath>
#include <tuple>
#include <unordered_set>

#include "types.h"
#include "bramble.h"
#include "g2t.h"
#include "evaluate.h"
#include "bam.h"
#include "mates.h"

#ifndef NOTHREADS
#include "GThreads.h"
#endif

#include "quill/Backend.h"
#include "quill/Frontend.h"
#include "quill/LogMacros.h"
#include "quill/Logger.h"
#include "quill/sinks/ConsoleSink.h"

#ifndef NOTHREADS
extern GFastMutex bam_io_mutex;     // protects BAM io
#endif
const uint32_t CHUNK_SIZE = 5000;   // size of each write to BAM
extern double similarity_threshold;

extern bool LONG_READS;
extern bool USE_FASTA;
extern bool SOFT_CLIPS;
extern bool QUIET;
extern bool BRAMBLE_DEBUG;

extern uint32_t total_complete;
extern uint32_t total_unique;
extern uint32_t dropped_reads;
extern uint32_t total_processed;
extern uint32_t print_mod;

extern uint32_t seen_last_out;

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
                            quill::Logger *logger,
                            uint8_t *seq, int seq_len) {

    // Evaluate read group to find matches
    std::vector<ExonChainMatch> matches = 
      evaluator->evaluate(read, g2t, seq, seq_len);

    total_processed++;
    if ((!QUIET || BRAMBLE_DEBUG) 
      && (total_processed % print_mod == 0)) {
      
      uint32_t out = int(total_processed / 100) * 100;
      if (out != seen_last_out) {
        LOG_INFO(logger, "{} alignments processed", out);
        seen_last_out = out;
      }
    }

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
        // genome CIGAR, kept for NM (update_cigar overwrites it); one buffer per thread
        thread_local std::vector<uint32_t> genome_cigar;
        genome_cigar.assign(new_cigar, new_cigar + n_cigar);
        int32_t nm = 0;
        update_cigar(b, new_cigar, n_cigar, cigar_mem, cigar, nm);
        // NM for the projected alignment (MD removed); before reverse-complementing
        set_projected_nm(b, genome_cigar, is_first ? pair->r_align : pair->m_align);
        // drop aux tags that describe the genome alignment (keep list in bam.cpp)
        keep_only_known_tags(b);

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
        order_tags(b);
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

  // 64-bit FNV-1a: a fixed, fully specified hash, so tie-breaking is the same
  // on every platform and standard library, and in the Rust version
  // (bramble-rs evaluate.rs deterministic_index)
  static uint64_t fnv1a_64(const std::string &s) {
    uint64_t h = 0xcbf29ce484222325ULL;
    for (unsigned char c : s) {
      h ^= c;
      h *= 0x100000001b3ULL;
    }
    return h;
  }

  // Sort key for tied candidates: transcript names and positions, not tids,
  // because tid numbering differs between the C++ and Rust versions; then the
  // input order of the genomic alignments (same in both versions), for ties on
  // the same transcript and position
  static std::tuple<std::string, pos_t, bool, std::string, pos_t, read_id_t, read_id_t>
  tie_key(const BamInfo *info, std::shared_ptr<g2tTree> &g2t) {
    auto pos = [](const AlignInfo &a) { return a.strand == '+' ? a.fwpos : a.rcpos; };
    read_id_t r_idx = info->read1 ? info->read1->index : 0;
    if (!info->is_paired) {
      return {g2t->getTidName(info->r_tid), pos(info->r_align), false, std::string(), 0,
              r_idx, 0};
    }
    read_id_t m_idx = info->read2 ? info->read2->index : 0;
    return {g2t->getTidName(info->r_tid), pos(info->r_align), true,
            g2t->getTidName(info->m_tid), pos(info->m_align), r_idx, m_idx};
  }

  // Output order of a read's units, the same in the Rust version (groups.rs
  // assign_primary_and_hit_indices): pairs before single records, then input
  // order of the genomic alignments, then transcript names and positions
  // (not tids: tid numbering differs between the two versions)
  static std::tuple<bool, read_id_t, read_id_t, std::string, pos_t, std::string, pos_t>
  unit_key(const BamInfo *info, std::shared_ptr<g2tTree> &g2t) {
    auto pos = [](const AlignInfo &a) { return a.strand == '+' ? a.fwpos : a.rcpos; };
    read_id_t r_idx = info->read1 ? info->read1->index : 0;
    if (!info->is_paired) {
      return {true, r_idx, 0, g2t->getTidName(info->r_tid), pos(info->r_align),
              std::string(), 0};
    }
    read_id_t m_idx = info->read2 ? info->read2->index : 0;
    return {false, r_idx, m_idx, g2t->getTidName(info->r_tid), pos(info->r_align),
            g2t->getTidName(info->m_tid), pos(info->m_align)};
  }

  // Read end of a record: 0 single-end, 1 read 1, 2 read 2
  static int end_of(const std::shared_ptr<ReadOut> &r) {
    uint16_t f = r->brec->get_b()->core.flag;
    if (!(f & BAM_FPAIRED)) return 0;
    return (f & BAM_FREAD1) ? 1 : ((f & BAM_FREAD2) ? 2 : 0);
  }

  // Identity of a projected record: its genomic alignment (read end, strand,
  // reference, position, CIGAR) and its transcript placement
  static std::string record_key(const std::shared_ptr<ReadOut> &r, tid_t tid,
                                const AlignInfo &a) {
    bam1_t* b = r->brec->get_b();
    uint16_t flag = b->core.flag & (BAM_FREAD1 | BAM_FREAD2 | BAM_FREVERSE);
    uint32_t n_cigar = b->core.n_cigar;
    pos_t pos = (a.strand == '+') ? a.fwpos : a.rcpos;
    std::string k;
    k.append((const char*)&flag, sizeof flag);
    k.append((const char*)&b->core.tid, sizeof b->core.tid);
    k.append((const char*)&b->core.pos, sizeof b->core.pos);
    k.append((const char*)&n_cigar, sizeof n_cigar);
    k.append((const char*)bam_get_cigar(b), n_cigar * sizeof(uint32_t));
    k.append((const char*)&tid, sizeof tid);
    k.append((const char*)&pos, sizeof pos);
    k.push_back(a.strand);
    return k;
  }

  // Collapse duplicate output units of one read: the same genomic alignment can
  // be reported once per mate (e.g. three identical read1 alignments with three
  // different mates), giving the same transcript record several times. Pairs
  // are kept first, then single records not already kept, each in input order,
  // so the earliest copy survives; dropped units are deleted. Mirrors Rust
  // collapse_duplicate_units (bramble-rs groups.rs). Call before NH/HI/primary.
  static void collapse_duplicates(std::vector<BamInfo*> &pairs) {
    auto idx = [](const std::shared_ptr<ReadOut> &r) { return r ? r->index : 0; };
    std::vector<BamInfo*> units(pairs);
    std::stable_sort(units.begin(), units.end(), [&](BamInfo* a, BamInfo* b) {
      return std::make_tuple(!a->is_paired, idx(a->read1),
                             a->is_paired ? idx(a->read2) : 0, a->r_tid)
           < std::make_tuple(!b->is_paired, idx(b->read1),
                             b->is_paired ? idx(b->read2) : 0, b->r_tid);
    });

    std::unordered_set<std::string> kept_records, kept_pairs;
    pairs.clear();
    for (BamInfo* u : units) {
      std::string a = record_key(u->read1, u->r_tid, u->r_align);
      bool keep;
      if (u->is_paired) {
        std::string b = record_key(u->read2, u->m_tid, u->m_align);
        keep = kept_pairs.insert(std::to_string(a.size()) + ':' + a + b).second;
        if (keep) { kept_records.insert(a); kept_records.insert(b); }
      } else {
        keep = kept_records.insert(a).second;
      }
      if (keep) pairs.push_back(u);
      else delete u;
    }
  }

  void convert_reads(std::vector<CReadAln> &reads,
                    std::shared_ptr<g2tTree> g2t, 
                    std::shared_ptr<ReadEvaluator> evaluator, 
                    quill::Logger *logger,
                    BamIO *io) {

    // Chunk buffer for BAM output
    std::vector<BamInfo*> bam_info;
    bam_info.reserve(CHUNK_SIZE*1.2);

    // Buffer for grouping by read name before filtering.
    // Must use std::string keys — const char* would use pointer comparison,
    // causing each genomic alignment of the same read to land in a separate
    // bucket and receive an incorrect per-alignment NH count.
    std::unordered_map<std::string, std::vector<BamInfo*>> pairs_by_name;
    std::vector<std::string> name_order;  // reads are written in input order
    uint32_t n_pairs = 0;

    auto flush = [&]() {
      // Filter by strand before writing
      std::vector<BamInfo*> filtered_bam_info;
      
      for (const std::string& read_name : name_order) {
        auto& pairs = pairs_by_name[read_name];
        // drop duplicate records (same genomic alignment, same placement)
        collapse_duplicates(pairs);
        // fixed unit order: HI follows it, and records are written in it
        std::stable_sort(pairs.begin(), pairs.end(), [&](BamInfo* a, BamInfo* b) {
          return unit_key(a, g2t) < unit_key(b, g2t);
        });

        // Determine primary/secondary alignment status
        auto best_it = pairs.end();
        double best_score = -std::numeric_limits<double>::infinity();
        int count_at_best = 0;
        int32_t hit_index = 1;
        int total_matches = 0;

        // HI per unit, shared by the two mates of a pair (as STAR); NH per read
        // end: the number of records of that end (SAM spec)
        uint32_t nh_by_end[3] = {0, 0, 0};
        for (auto it = pairs.begin(); it != pairs.end(); ++it) {
          auto &info = (*it);
          info->r_align.hit_index = hit_index;
          nh_by_end[end_of(info->read1)]++;
          total_matches++;

          if (info->is_paired) {
            info->m_align.hit_index = hit_index;
            nh_by_end[end_of(info->read2)]++;
            total_matches++;
          }
          hit_index++;

          // Best score for this pair is the higher of the two
          double pair_score = info->r_align.similarity_score;
          if (info->is_paired) {
            pair_score = std::max(pair_score, info->m_align.similarity_score);
          }

          if (pair_score > best_score) {
            best_score = pair_score;
            best_it = it;
            count_at_best = 1;
          } else if (pair_score == best_score) {
            count_at_best++;
          }
        }

        if (best_it != pairs.end()) {
          if (count_at_best == 1) {
            auto &best_info = (*best_it);
            best_info->r_align.primary_alignment = true;
            if (best_info->is_paired) {
              best_info->m_align.primary_alignment = true;
            }
          } else {
            // Tie: pick deterministically among the tied candidates, sorted by
            // transcript name and position, using a fixed hash of the read name,
            // so the choice is the same on every platform and in the Rust version
            std::vector<decltype(pairs.begin())> tied_its;
            for (auto it = pairs.begin(); it != pairs.end(); ++it) {
              auto &info = (*it);
              double pair_score = info->r_align.similarity_score;
              if (info->is_paired) {
                pair_score = std::max(pair_score, info->m_align.similarity_score);
              }
              if (pair_score == best_score) {
                tied_its.push_back(it);
              }
            }

            std::sort(tied_its.begin(), tied_its.end(), [&](auto a, auto b) {
              return tie_key(*a, g2t) < tie_key(*b, g2t);
            });
            size_t idx = fnv1a_64(read_name) % tied_its.size();
            auto &secondary = *(tied_its[idx]);

            secondary->r_align.primary_alignment = true;
            if (secondary->is_paired) {
              secondary->m_align.primary_alignment = true;
            }
          }
        }

        // Every read end gets one primary. If the best unit is a single record,
        // the other end's best record (any unit) becomes primary on its own; its
        // partner, if any, stays secondary. Ties: sorted by transcript name,
        // position and input order, picked by the read-name hash (as Rust
        // assign_primary_and_hit_indices).
        {
          struct Rec { AlignInfo* align; tid_t tid; read_id_t index; int end; };
          std::vector<Rec> recs;
          std::vector<int> covered;
          for (auto* info : pairs) {
            recs.push_back({&info->r_align, info->r_tid, info->read1->index, end_of(info->read1)});
            if (info->is_paired)
              recs.push_back({&info->m_align, info->m_tid, info->read2->index, end_of(info->read2)});
          }
          for (const auto &r : recs)
            if (r.align->primary_alignment) covered.push_back(r.end);
          for (int end : {1, 2}) {
            if (std::find(covered.begin(), covered.end(), end) != covered.end()) continue;
            double best = -std::numeric_limits<double>::infinity();
            std::vector<const Rec*> cands;
            for (const auto &r : recs) {
              if (r.end != end) continue;
              if (r.align->similarity_score > best) { best = r.align->similarity_score; cands.clear(); }
              if (r.align->similarity_score == best) cands.push_back(&r);
            }
            if (cands.empty()) continue;
            auto key = [&](const Rec* r) {
              pos_t p = (r->align->strand == '+') ? r->align->fwpos : r->align->rcpos;
              return std::make_tuple(g2t->getTidName(r->tid), p, r->index);
            };
            std::sort(cands.begin(), cands.end(),
                      [&](const Rec* a, const Rec* b) { return key(a) < key(b); });
            cands[fnv1a_64(read_name) % cands.size()]->align->primary_alignment = true;
          }
        }

        // Update statistics
        total_complete += total_matches;
        if (pairs.size() == 1) total_unique++;

        // NH and MAPQ per read end
        for (auto* pair : pairs) {
          pair->read1->nh = nh_by_end[end_of(pair->read1)];
          pair->read1->mapq = get_mapq(pair->read1->nh);
          if (pair->is_paired && pair->read2) {
            pair->read2->nh = nh_by_end[end_of(pair->read2)];
            pair->read2->mapq = get_mapq(pair->read2->nh);
          }
          filtered_bam_info.push_back(pair);
        }
      }
      
      write_to_bam(io, filtered_bam_info);
      for (auto* info : filtered_bam_info) delete info;
      pairs_by_name.clear();
      name_order.clear();
      n_pairs = 0;
    };

    auto emit_pair = [&](BamInfo* pair, bool /*is_last*/) {
      if (pair && pair->read1) {
        auto [it, inserted] = pairs_by_name.try_emplace(pair->read1->brec->name());
        if (inserted) name_order.push_back(it->first);
        it->second.push_back(pair);
        n_pairs++;
        // Flushing is deferred to the end of each read-name group so that all
        // genomic alignments of a multi-mapper are grouped before NH is computed.
      }
    };

    // Track seen reads
    unordered_set<read_id_t> seen;
    seen.reserve(reads.size() * 1.2);

    for (read_id_t id = 0; id < reads.size(); ) {

      // Group reads with identical read names
      read_id_t start = id;
      const std::string& name = reads[id].brec->name();

      int seq_len = 0;
      uint8_t* seq = nullptr;
      bool found = false;

      if (USE_FASTA) {
        bam1_t* b = reads[id].brec->get_b();
        seq_len = b->core.l_qseq;
        if (seq_len > 0) {      
          seq = bam_get_seq(b);
          found = true;
        }
      }

      id++;

      while (id < reads.size() && reads[id].brec->name() == name) {
        if (USE_FASTA && !found) {
          bam1_t* b = reads[id].brec->get_b();
          seq_len = b->core.l_qseq;
          if (seq_len > 0) { 
            seq = bam_get_seq(b);
            found = true;
          }
        }
        id++;
      }

      read_id_t end = id; // [start, end)
      
      // Process read groups
      bool dropped = true;
      for (read_id_t i = start; i < end; i++) {
        if (seen.count(i)) continue;

        ReadInfo* this_read = process_read_out(reads[i], 
          i, g2t, evaluator, logger, seq, seq_len);
        if (this_read) dropped = false;

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

          ReadInfo* mate_read = process_read_out(reads[mate_id], 
            mate_id, g2t, evaluator, logger, seq, seq_len);

          if (mate_read) dropped = false;
          process_mate_pair(this_read, mate_read, emit_pair);
          delete mate_read;
          seen.insert(mate_id);
        }

        delete this_read;
        seen.insert(i);
      }

      if (dropped) dropped_reads++;

      // All genomic alignments for this read name have been emitted.
      // Flush the accumulated buffer if it has grown large enough.
      if (n_pairs >= CHUNK_SIZE) {
        flush();
      }
    }

    if (!pairs_by_name.empty()) flush();
  }

}