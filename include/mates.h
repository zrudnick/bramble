
#pragma once

#include "types.h"

using namespace bramble;

namespace bramble {

  struct ReadInfo;
  struct BamInfo;
  struct BundleData;

  void add_mate_info(const std::unordered_set<tid_t> &final_transcripts,
                    const std::unordered_set<tid_t> &read_transcripts,
                    const std::unordered_set<tid_t> &mate_transcripts,
                    const std::unordered_map<tid_t, pos_t> &read_positions,
                    const std::unordered_map<tid_t, pos_t> &mate_positions,
                    ReadInfo* this_read, ReadInfo* mate_read, uint8_t mate_case, 
                    std::function<void(BamInfo*, bool)> emit_pair);
                    
  void update_read_matches(ReadInfo *read_info,
                          const std::unordered_set<tid_t> &final_transcripts);

  void process_mate_pair(ReadInfo* this_read, ReadInfo* mate_read,
                        std::function<void(BamInfo*, bool)> emit_pair);
}