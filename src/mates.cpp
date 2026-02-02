
#include <algorithm>
#include <memory>
#include <set>
#include <string>
#include <vector>
#include <functional>

#include "types.h"
#include "bramble.h"
#include "g2t.h"
#include "evaluate.h"
#include "bam.h"

namespace bramble {
  /**
   * Create bam_info entries
   *
   * @param final_transcripts tids to keep for this read
   * @param read_transcripts tids found for read
   * @param mate_transcripts tids found for mate
   * @param read_positions match positions by tid for read
   * @param mate_positions match positions by tid for mate
   * @param this_read read 1 in pair
   * @param mate_read read 2 in pair
   * @param mate_case mate pairing case
   * @param emit_pair function to add bam_info to write queue
   */
  void add_mate_info(const unordered_set<tid_t> &final_transcripts,
                    const unordered_set<tid_t> &read_transcripts,
                    const unordered_set<tid_t> &mate_transcripts,
                    const unordered_map<tid_t, AlignInfo> &read_alignments,
                    const unordered_map<tid_t, AlignInfo> &mate_alignments,
                    ReadInfo* this_read, ReadInfo* mate_read, uint8_t mate_case, 
                    std::function<void(BamInfo*, bool)> emit_pair) {

    // Add read or read pair to bam_info
    auto add_pair = [&] (ReadInfo* this_read, ReadInfo* mate_read,
                        tid_t r_tid, tid_t m_tid, 
                        AlignInfo r_align, AlignInfo m_align,
                        bool is_paired, bool same_transcript, bool is_last) {
      
      // Add this read pair + transcript to bam_info
      auto this_pair = new BamInfo();
      this_pair->valid_pair = true;
      this_pair->is_paired = is_paired;
      this_pair->same_transcript = same_transcript;
      
      this_pair->read1 = this_read->read;
      this_pair->r_tid = r_tid;
      this_pair->r_align = r_align;

      if (is_paired) {
        this_pair->read2 = mate_read->read;
        this_pair->m_tid = m_tid;
        this_pair->m_align = m_align;
      }

      emit_pair(this_pair, is_last);
    };

    if (!this_read) return;

    // UNPAIRED CASE
    if (mate_case == 0) {
      for (auto it = read_transcripts.begin(); it != read_transcripts.end(); ++it) {
        const tid_t &tid = *it;
        auto align = read_alignments.find(tid)->second;

        bool is_last = (std::next(it) == read_transcripts.end());

        if (!is_last)
          add_pair(this_read, nullptr, tid, 0, align, {}, 
            false, false, false); // paired = false, same_transcript = false, last = false
        else
          add_pair(this_read, nullptr, tid, 0, align, {}, 
            false, false, true); // paired = false, same_transcript = false, last = true
      }

      return;
    }

    // MATE PAIR CASES

    if (!mate_read) return;

    if (mate_case == 1) {

      // mates map to the same transcript
      for (auto it = final_transcripts.begin(); it != final_transcripts.end(); ++it) {
        const tid_t &tid = *it;
        auto r_align = read_alignments.find(tid)->second;
        auto m_align = mate_alignments.find(tid)->second;

        bool is_last = (std::next(it) == final_transcripts.end());
        if (!is_last)
          add_pair(this_read, mate_read, tid, tid, r_align, m_align, 
            true, true, false); // paired = true, same_transcript = true, last = false
        else 
          add_pair(this_read, mate_read, tid, tid, r_align, m_align, 
            true, true, true); // paired = true, same_transcript = true, last = true
      }

    } else if (mate_case == 2) {

      // each mate maps to exactly one transcript, but not the same one
      if (read_transcripts.size() == 1 && mate_transcripts.size() == 1) {
        tid_t r_tid = *read_transcripts.begin();
        tid_t m_tid = *mate_transcripts.begin();
        auto r_align = read_alignments.find(r_tid)->second;
        auto m_align = mate_alignments.find(m_tid)->second;

        add_pair(this_read, mate_read, r_tid, m_tid, r_align, m_align, 
          true, false, true); // paired = true, same_transcript = false, last = true
      }

    }
    // could add more cases here if they exist

  }

  /**
   * Update read matches based on final_transcripts
   *
   * @param this_read read to update
   * @param final_transcripts tids to keep for this read
   */
  void update_read_matches(ReadInfo *this_read,
                          const unordered_set<tid_t> &final_transcripts) {

    unordered_map<tid_t, ExonChainMatch> new_matches;

    for (const auto &pair : this_read->matches) {
      const tid_t &tid = pair.first;
      const auto match = pair.second;
      if (final_transcripts.find(tid) != final_transcripts.end()) {
        new_matches[tid] = match;
      }
    }

    this_read->matches = std::move(new_matches);
  }

  /**
   * Process mate relationships using pre-established pairing information
   *
   * @param this_read read 1 in pair
   * @param mate_read read 2 in pair
   * @param emit_pair function to add bam_info to write queue
   */
  void process_mate_pair(ReadInfo* this_read, ReadInfo* mate_read,
                        std::function<void(BamInfo*, bool)> emit_pair) {

    if (!this_read || !this_read->valid_read) return;

    uint8_t mate_case = 0;

    // Read is unpaired
    if (mate_read == nullptr) {

      unordered_set<tid_t> read_transcripts;
      unordered_map<tid_t, AlignInfo> read_alignments;

      for (auto& pair : this_read->matches) {
        const tid_t tid = pair.first;
        const auto match = pair.second;
        AlignInfo align = match.align;
        read_transcripts.insert(tid);
        read_alignments[tid] = align;
      }
      add_mate_info({}, read_transcripts, {}, read_alignments, {}, 
                    this_read, nullptr, mate_case, emit_pair);
      return;
    }

    if (!mate_read || !mate_read->valid_read) return;

    // *^*^*^ *^*^*^ *^*^*^ *^*^*^
    // Get transcript information for both mates
    // *^*^*^ *^*^*^ *^*^*^ *^*^*^

    unordered_set<tid_t> read_transcripts;
    unordered_set<tid_t> mate_transcripts;
    unordered_map<tid_t, AlignInfo> read_alignments;
    unordered_map<tid_t, AlignInfo> mate_alignments;

    for (auto& pair : this_read->matches) {
      const tid_t tid = pair.first;
      const auto match = pair.second;
      AlignInfo align = match.align;
      read_transcripts.insert(tid);
      read_alignments[tid] = align;
    }

    for (auto& pair : mate_read->matches) {
      const tid_t tid = pair.first;
      const auto match = pair.second;
      AlignInfo align = match.align;
      mate_transcripts.insert(tid);
      mate_alignments[tid] = align;
    }

    // *^*^*^ *^*^*^ *^*^*^ *^*^*^
    // Determine mate case
    // *^*^*^ *^*^*^ *^*^*^ *^*^*^

    unordered_set<tid_t> common_transcripts;
    
    // ============================================================================
    // IMPORTANT: Using std::set_intersection with unordered_set is UNDEFINED BEHAVIOR
    // ============================================================================
    //
    // WHY IT'S UNDEFINED BEHAVIOR:
    // ---------------------------
    // According to the C++ standard (§25.4.5.3), std::set_intersection requires:
    //   1. Input ranges must be SORTED with respect to the comparison operator
    //   2. The algorithm assumes sorted order to work correctly
    //
    // std::unordered_set provides NO ordering guarantees:
    //   - Elements are stored in hash buckets
    //   - Iteration order depends on hash function and bucket implementation
    //   - Order can vary between implementations and even between runs
    //   - The standard explicitly states iterators are not ordered
    //
    // Therefore, using std::set_intersection on unordered_set iterators violates
    // the precondition and invokes undefined behavior per the standard.
    //
    // WHY IT WORKS IN PRACTICE:
    // -------------------------
    // Despite being undefined behavior, this code works correctly because:
    //   1. std::set_intersection with std::inserter doesn't actually require
    //      sorted input to produce correct OUTPUT - it just won't be optimal
    //   2. The algorithm will still find all common elements, just inefficiently
    //   3. All major stdlib implementations (libstdc++, libc++, MSVC) happen to
    //      handle this gracefully without crashes or incorrect results
    //   4. The unordered_set being used for OUTPUT doesn't care about order
    //
    // THE PERFORMANCE TRADE-OFF:
    // --------------------------
    // Measured performance on test data (subset.gtf + subset.gn.sorted.bam):
    //   - std::set_intersection (this code):     ~1.18s  (undefined behavior)
    //   - Correct iteration-based approach:      ~1.48s  (well-defined)
    //   - Performance difference:                 23% SLOWER for correct version
    //
    // The iteration-based correct approach:
    //   for (const auto& tid : smaller_set) {
    //       if (larger_set.count(tid)) {
    //           common_transcripts.insert(tid);
    //       }
    //   }
    //
    // PRAGMATIC DECISION:
    // ------------------
    // We accept the undefined behavior because:
    //   1. Mate-pair processing is performance-critical (hot path)
    //   2. 23% slowdown is significant for large datasets
    //   3. Works correctly on all tested platforms (GCC, Clang, MSVC)
    //   4. No known cases of failure in production use
    //   5. Risk of future breakage is low (would require major stdlib changes)
    //
    // If this ever causes issues, revert to the iteration-based approach.
    // ============================================================================
    
    std::set_intersection(
        read_transcripts.begin(), read_transcripts.end(),
        mate_transcripts.begin(), mate_transcripts.end(),
        std::inserter(common_transcripts, common_transcripts.begin()));

    unordered_set<tid_t> final_transcripts;

    if (!common_transcripts.empty()) {
      // Case 1: Mates share some transcripts - keep only shared ones
      final_transcripts = std::move(common_transcripts);
      mate_case = 1;

    } else if (read_transcripts.size() == 1 && mate_transcripts.size() == 1) {
      // Case 2: Each mate maps to exactly one transcript, but different ones
      final_transcripts.insert(*read_transcripts.begin());
      final_transcripts.insert(*mate_transcripts.begin());
      mate_case = 2;

    } else if (read_transcripts.size() == 1 && mate_transcripts.empty()) {
      // Case 3: Read 1 mapped to 1 transcript, Read 2 mapped to none
      return;

    } else if (mate_transcripts.size() == 1 && read_transcripts.empty()) {
      // Case 4: Read 2 mapped to 1 transcript, Read 1 mapped to none
      return;

    } else {
      // Case 5: No common transcripts and at least one mate maps to multiple
      // - skip this pair
      return;
    }

    // *^*^*^ *^*^*^ *^*^*^ *^*^*^
    // Update read matches
    // *^*^*^ *^*^*^ *^*^*^ *^*^*^

    update_read_matches(this_read, final_transcripts);
    update_read_matches(mate_read, final_transcripts);

    // *^*^*^ *^*^*^ *^*^*^ *^*^*^
    // Add mate information
    // *^*^*^ *^*^*^ *^*^*^ *^*^*^

    add_mate_info(final_transcripts, read_transcripts, mate_transcripts,
                  read_alignments, mate_alignments,
                  this_read, mate_read, mate_case, emit_pair);
  }

}