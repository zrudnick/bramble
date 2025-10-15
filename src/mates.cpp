
#include <algorithm>
#include <memory>
#include <set>
#include <string>
#include <vector>
#include <unordered_set>
#include <functional>

#include "types.h"
#include "bundles.h"
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
  void add_mate_info(const std::unordered_set<tid_t> &final_transcripts,
                    const std::unordered_set<tid_t> &read_transcripts,
                    const std::unordered_set<tid_t> &mate_transcripts,
                    const std::unordered_map<tid_t, pos_t> &read_positions,
                    const std::unordered_map<tid_t, pos_t> &mate_positions,
                    ReadInfo* this_read, ReadInfo* mate_read, uint8_t mate_case, 
                    std::function<void(BamInfo*, bool)> emit_pair) {

    // Add read or read pair to bam_info
    auto add_pair = [&] (ReadInfo* this_read, ReadInfo* mate_read,
                        tid_t r_tid, tid_t m_tid, 
                        pos_t r_pos, pos_t m_pos, 
                        bool is_paired, bool is_last) {
      
      // Add this read pair + transcript to bam_info
      auto this_pair = new BamInfo();
      this_pair->valid_pair = true;
      this_pair->is_paired = is_paired;
      
      this_pair->r_tid = r_tid;
      this_pair->r_pos = r_pos;
      this_pair->read1 = this_read->read;

      if (is_paired) {
        this_pair->read2 = mate_read->read;
        this_pair->m_tid = m_tid;
        this_pair->m_pos = m_pos;
      }

      emit_pair(this_pair, is_last);
    };

    if (!this_read) return;

    // UNPAIRED CASE
    if (mate_case == 0) {
      for (auto it = read_transcripts.begin(); it != read_transcripts.end(); ++it) {
        const tid_t &tid = *it;
        auto r_it = read_positions.find(tid);
        auto pos = r_it->second;

        bool is_last = (std::next(it) == final_transcripts.end());

        if (!is_last)
          add_pair(this_read, nullptr, tid, 0, pos, 0, false, false); 
          // paired = false, last = false
        else
          add_pair(this_read, nullptr, tid, 0, pos, 0, false, true); 
          // paired = false, last = true
      }

      return;
    }

    // MATE PAIR CASES

    if (!mate_read) return;

    if (mate_case == 1) {
      for (auto it = final_transcripts.begin(); it != final_transcripts.end(); ++it) {
        const tid_t &tid = *it;
        auto r_it = read_positions.find(tid);
        auto m_it = mate_positions.find(tid);
        pos_t r_pos = r_it->second;
        pos_t m_pos = m_it->second;

        bool is_last = (std::next(it) == final_transcripts.end());
        if (!is_last)
          add_pair(this_read, mate_read, tid, tid, r_pos, m_pos, true, false); 
          // paired = true, last = false
        else 
          add_pair(this_read, mate_read, tid, tid, r_pos, m_pos, true, true); 
          // paired = true, last = true
      }

    } else if (mate_case == 2) {

      // each mate maps to exactly one transcript, but not the same one
      if (read_transcripts.size() == 1 && mate_transcripts.size() == 1) {
        tid_t r_tid = *read_transcripts.begin();
        tid_t m_tid = *mate_transcripts.begin();
        pos_t r_pos = read_positions.find(r_tid)->second;
        pos_t m_pos = mate_positions.find(m_tid)->second;

        add_pair(this_read, mate_read, r_tid, m_tid, r_pos, m_pos, true, true); 
        // paired = true, last = true
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
                          const std::unordered_set<tid_t> &final_transcripts) {

    std::vector<ExonChainMatch> new_matches;

    for (const auto &match : this_read->matches) {
      const tid_t &tid = match.tid;
      if (final_transcripts.find(tid) != final_transcripts.end()) {
        new_matches.emplace_back(match);
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

      std::unordered_set<tid_t> read_transcripts;
      std::unordered_map<tid_t, pos_t> read_positions;

      for (auto& match : this_read->matches) {
        tid_t tid = match.tid;
        uint pos = match.pos;
        read_transcripts.insert(tid);
        read_positions[tid] = pos;
      }
      add_mate_info({}, read_transcripts, {}, read_positions, {}, 
                    this_read, nullptr, mate_case, emit_pair);
      return;
    }

    if (!mate_read || !mate_read->valid_read) return;

    // *^*^*^ *^*^*^ *^*^*^ *^*^*^
    // Get transcript information for both mates
    // *^*^*^ *^*^*^ *^*^*^ *^*^*^

    std::unordered_set<tid_t> read_transcripts;
    std::unordered_set<tid_t> mate_transcripts;
    std::unordered_map<tid_t, pos_t> read_positions;
    std::unordered_map<tid_t, pos_t> mate_positions;

    for (auto& match : this_read->matches) {
      tid_t tid = match.tid;
      uint pos = match.pos;
      read_transcripts.insert(tid);
      read_positions[tid] = pos;
    }

    for (auto& match : mate_read->matches) {
      tid_t tid = match.tid;
      uint pos = match.pos;
      mate_transcripts.insert(tid);
      mate_positions[tid] = pos;
    }

    // *^*^*^ *^*^*^ *^*^*^ *^*^*^
    // Determine mate case
    // *^*^*^ *^*^*^ *^*^*^ *^*^*^

    std::unordered_set<tid_t> common_transcripts;
    std::set_intersection(
        read_transcripts.begin(), read_transcripts.end(),
        mate_transcripts.begin(), mate_transcripts.end(),
        std::inserter(common_transcripts, common_transcripts.begin()));

    std::unordered_set<tid_t> final_transcripts;

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
                  read_positions, mate_positions,
                  this_read, mate_read, mate_case, emit_pair);
  }

}