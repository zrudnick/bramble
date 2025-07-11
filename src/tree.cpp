
#include <algorithm>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "bam.h"
#include "bramble.h"

extern bool USE_FASTA;
uint32_t overhang_threshold = 8;
const uint32_t max_mappings_per_read = 1000;
using tid_t = uint32_t;
using read_id_t = uint32_t; // read index in bundle
using bam_id_t = uint32_t;  // id in bam_info

/**
 * Prints g2t tree with nodes sorted by genome coordinate
 * For debugging only
 *
 * @param g2t g2t tree for bundle
 */
void print_tree(g2tTree *g2t) {
  auto intervals = g2t->fw_tree->getOrderedIntervals();
  printf("FW Tree:\n");
  for (auto &interval : intervals) {
    printf("(%u, %u)", interval->start, interval->end);
    for (auto &tid : interval->tids) {
      printf(" %d ", tid);
      printf(" %d ", interval->tid_cum_len[tid]);
    }
    printf("-> ");
  }
  printf("END\n");

  intervals = g2t->rc_tree->getOrderedIntervals();
  printf("RC Tree:\n");
  for (auto &interval : intervals) {
    printf("(%u, %u)", interval->start, interval->end);
    for (auto &tid : interval->tids) {
      printf(" %d ", tid);
      printf(" %d ", interval->tid_cum_len[tid]);
    }
    printf("-> ");
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
    tid_t tid = g2t->insertTidString(tid_string, io);

    for (int j = 0; j < guide->exons.Count(); j++) {
      uint g_start, g_end;
      if (strand == '+') {
        g_start = guide->exons[j]->start - bundle->start;
        g_end = guide->exons[j]->end - bundle->start;
      } else {
        g_start = bundle->end - guide->exons[j]->end;
        g_end = bundle->end - guide->exons[j]->start;
      }

      g2t->addGuideExon(g_start, g_end, tid, strand);
    }
  }

  g2t->buildAllTidChains();
  g2t->precomputeAllCumulativeLengths();

  // print_tree(g2t);
  return g2t;
}

std::string reverse_complement(const std::string &seq) {
  std::string rc;
  rc.reserve(seq.length());

  for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
    switch (*it) {
    case 'A':
      rc += 'T';
      break;
    case 'T':
      rc += 'A';
      break;
    case 'G':
      rc += 'C';
      break;
    case 'C':
      rc += 'G';
      break;
    case 'N':
      rc += 'N';
      break;
    default:
      rc += *it;
      break;
    }
  }
  return rc;
}

std::string extract_sequence(const char *gseq, uint start, uint length,
                             char strand) {
  if (gseq == nullptr)
    return "";

  // Get substring (or reverse complement)
  std::string seq(gseq + start, length);
  if (strand == '-') {
    seq = reverse_complement(seq);
  }
  return seq;
}

// Function to check backward overhang (left extension)
bool check_backward_overhang(IntervalNode *interval, uint exon_start,
                             char read_strand,
                             std::set<tid_t> &exon_tids, g2tTree *g2t,
                             BundleData *bundle) {

  std::set<tid_t> tmp_exon_tids;
  uint overhang_start;
  uint overhang_length = interval->start - exon_start;
  if (overhang_length > overhang_threshold)
    return false; // overhang is too long

  if (read_strand == '+') {
    overhang_start = exon_start;
  } else {
    uint relative_end = bundle->end - bundle->start;
    overhang_start = relative_end - exon_start - overhang_length;
  }

  // Extract overhang sequence from genomic sequence (left extension)
  std::string overhang_seq = extract_sequence(bundle->gseq, overhang_start,
                                              overhang_length, read_strand);
  for (const auto &tid : exon_tids) {

    // Use the TID-based edge to get the previous interval directly from current
    // interval
    IntervalNode *prev_interval = g2t->getPrevNode(interval, tid, read_strand);

    if (prev_interval &&
        prev_interval->end <= overhang_start + overhang_length) {
      // This is the previous guide exon for our transcript
      // Extract sequence from the end of this guide exon
      uint check_length =
          std::min(overhang_length, prev_interval->end - prev_interval->start);
      uint prev_overhang_start;

      if (read_strand == '+') {
        prev_overhang_start = prev_interval->end - check_length;
      } else {
        prev_overhang_start = bundle->end - prev_interval->end - bundle->start;
      }

      std::string guide_seq = extract_sequence(
          bundle->gseq, prev_overhang_start, check_length, read_strand);

      if (overhang_seq == guide_seq) {
        tmp_exon_tids.insert(tid);
      }
    }
  }
  std::swap(exon_tids, tmp_exon_tids);
  tmp_exon_tids.clear();

  // Return true if we found matching transcripts
  return !exon_tids.empty();
}

std::set<tid_t> collapse_intervals(std::vector<IntervalNode *> sorted_intervals,
                                   uint exon_start, 
                                   bool is_first_exon, 
                                   char read_strand, 
                                   g2tTree *g2t,
                                   BundleData *bundle,
                                   bool &used_backwards_overhang,
                                   uint &soft_clip,
                                   IntervalNode *prev_last_interval
                                   ) {

  // 1. Ensure each read exon is fully covered by a series of adjacent guide
  // exons
  std::set<tid_t> exon_tids;

  auto first_interval = sorted_intervals[0];
  uint prev_end = first_interval->end;

  // PROCESS FIRST INTERVAL: LATER EXONS

  if (!is_first_exon) {
    // Copy TIDs from interval
    // --> skip any with exons between current and previous interval
    for (const auto &tid : first_interval->tids) {
      auto tid_last_interval =
          g2t->getPrevNode(first_interval, tid, read_strand);
      if ((tid_last_interval == prev_last_interval) 
          && (first_interval->start <= exon_start)){
        exon_tids.insert(tid);
      }
    }

  // PROCESS FIRST INTERVAL: FIRST EXON

    // Start is within bounds --> copy all TIDs from interval
  } else if (first_interval->start <= exon_start) {
    exon_tids.insert(first_interval->tids.begin(), first_interval->tids.end());

    // Use input FASTA (-S) to check for previous transcript exons to map to
  } else if (USE_FASTA) {
    // Copy TIDs
    exon_tids.insert(first_interval->tids.begin(), first_interval->tids.end());
    used_backwards_overhang = check_backward_overhang(
        first_interval, exon_start, read_strand, exon_tids, g2t, bundle);
    if (!used_backwards_overhang) {
      soft_clip = first_interval->start - exon_start;
    }

    // Soft clip backwards
  } else {
    soft_clip = first_interval->start - exon_start;
    return {};
  }

  // PROCESS FOLLOWING INTERVALS

  // Following intervals: check for TIDs, remove any we don't see
  for (size_t i = 1; i < sorted_intervals.size(); ++i) {
    const auto &interval = sorted_intervals[i];

    // Check continuity
    if (interval->start != prev_end)
      return {};

    // Remove TIDs not found in current interval
    auto it = exon_tids.begin();
    while (it != exon_tids.end()) {
      if (std::find(interval->tids.begin(), interval->tids.end(), *it) ==
          interval->tids.end()) {
        it = exon_tids.erase(it);
      } else {
        ++it;
      }
    }

    // Early exit if no common TIDs remain
    if (exon_tids.empty())
      return {};

    prev_end = interval->end;
  }

  return exon_tids;
}

/**
 * Check if we can match forwards to another interval (USE_FASTA mode)
 *
 * @param interval first interval the read matched to
 */
bool check_forward_overhang(IntervalNode *interval, uint exon_end,
                            char read_strand,
                            std::set<tid_t> &exon_tids, g2tTree *g2t,
                            BundleData *bundle) {

  std::set<tid_t> tmp_exon_tids;
  uint overhang_start;
  uint overhang_length = exon_end - interval->end;
  if (overhang_length > overhang_threshold)
    return false; // overhang is too long

  if (read_strand == '+') {
    overhang_start = interval->end;
  } else {
    uint relative_end = bundle->end - bundle->start;
    overhang_start = relative_end - interval->end - overhang_length;
  }

  // Extract overhang sequence from genomic sequence
  std::string overhang_seq = extract_sequence(bundle->gseq, overhang_start,
                                              overhang_length, read_strand);

  // For each transcript that was compatible so far
  for (const auto &tid : exon_tids) {

    // Find the next guide exon for this transcript after position i
    IntervalNode *next_interval = g2t->getNextNode(interval, tid, read_strand);

    if (next_interval && next_interval->start >= overhang_start) {
      // This is the next guide exon for our transcript
      // Extract sequence from the start of this guide exon
      uint check_length =
          std::min(overhang_length, next_interval->end - next_interval->start);
      uint next_overhang_start;

      if (read_strand == '+' || read_strand == 1) {
        next_overhang_start = next_interval->start;
      } else {
        next_overhang_start =
            bundle->end - next_interval->start - check_length - bundle->start;
      }

      std::string guide_seq = extract_sequence(
          bundle->gseq, next_overhang_start, check_length, read_strand);

      if (overhang_seq == guide_seq) {
        tmp_exon_tids.insert(tid);
      }
    }
  }
  
  std::swap(exon_tids, tmp_exon_tids);
  tmp_exon_tids.clear();

  // Return true if we found matching transcripts
  return !exon_tids.empty();
}

/**
 * Get the match position for this transcript
 *
 * @param interval first interval the read matched to
 * @param tid transcript id
 * @param read_strand strand that the read aligned to
 * @param g2t g2t tree for bundle
 * @param exon_start start of first read exon
 * @param used_backwards_overhang did we match backwards to a previous interval?
 * (USE_FASTA mode)
 */
uint get_match_pos(IntervalNode *interval, tid_t tid, char read_strand,
                   g2tTree *g2t, uint exon_start,
                   bool used_backwards_overhang) {
  IntervalNode *first_interval = interval;

  if (used_backwards_overhang) {
    first_interval = g2t->getPrevNode(interval, tid, read_strand);
    GMessage("used backwards overhang\n");
  }

  uint prev_node_sum =
      g2t->getCumulativeLength(first_interval, tid, read_strand);
  return (exon_start - first_interval->start) + prev_node_sum;
}

/**
 * Process a single read - returns valid matches
 *
 * @param read current read
 * @param bundle current bundle of reads
 * @param g2t g2t tree for bundle
 * @param bam_info read_info for all reads in bundle
 * @param read_index index of read in bundle->reads
 */
void process_read_out(BundleData *&bundle, 
                      uint read_index, 
                      g2tTree *g2t,
                      std::map<read_id_t, ReadInfo *> &read_info,
                      std::vector<uint32_t> group
                      ) {

  CReadAln *read = bundle->reads[read_index];
  std::set<std::tuple<tid_t, uint>> matches;
  std::set<std::tuple<tid_t, uint>> tmp_matches;

  GVec<GSeg> read_exons = read->segs;
  std::string read_name = read->brec->name();
  char read_strand = read->strand;

  IntervalNode *first_interval = nullptr;
  IntervalNode *last_interval = nullptr;
  IntervalNode *prev_last_interval = nullptr;

  uint exon_count = read_exons.Count();

  for (uint j = 0; j < exon_count; j++) {
    GSeg curr_exon = read_exons[j];
    uint exon_start, exon_end;

    std::vector<IntervalNode*> sorted_intervals;

    // TODO: Improve method of inferring splice strand
    if (read_strand == '+' || read_strand == '.') {
      if (read_strand == '.') read_strand = '+'; // this should only happen for first read exon
      exon_start = curr_exon.start - bundle->start;
      exon_end = curr_exon.end - bundle->start;
    } else {
      exon_start = bundle->end - curr_exon.end;
      exon_end = bundle->end - curr_exon.start;
    }

    // Find all guide intervals that contain the read exon
    sorted_intervals = g2t->getIntervals(exon_start, exon_end, read_strand);

     // If uncertain strand and forward didn't work, try reverse
    if (read_strand == '.' && sorted_intervals.empty()) { // this should only happen for first read exon
      read_strand = '-';
      sorted_intervals = g2t->getIntervals(exon_start, exon_end, read_strand);
    }
    
    // No overlap at all
    if (sorted_intervals.empty())
      return;

    // NOTE: This must come *after* the empty check!!
    prev_last_interval = last_interval;
    first_interval = sorted_intervals[0];
    last_interval = sorted_intervals[sorted_intervals.size() - 1];

    bool is_first_exon = (j == 0) ? true : false;
    bool is_last_exon = (j == exon_count - 1) ? true : false;

    bool used_backwards_overhang = false;
    uint soft_clip_front; // number bases to soft clip at front
    uint soft_clip_back;  // number bases to soft clip at back

    // *^*^*^ *^*^*^ *^*^*^ *^*^*^
    // First, check the intervals for this read exon
    // *^*^*^ *^*^*^ *^*^*^ *^*^*^

    auto exon_tids = collapse_intervals(
        sorted_intervals, exon_start, is_first_exon, read_strand, g2t, bundle,
        used_backwards_overhang, soft_clip_front, prev_last_interval);

    // Exon extends beyond guide interval (USE_FASTA mode)
    if (USE_FASTA && (is_last_exon) && (last_interval->end < exon_end)) {
      // Update valid transcripts to only those with matching extensions
      if (!check_forward_overhang(last_interval, exon_end, read_strand,
                                  exon_tids, g2t, bundle))
        soft_clip_back = exon_end - last_interval->end;

      // Last exon, Exon extends beyond guide interval (no fasta)
    } else if ((is_last_exon) && last_interval->end < exon_end) {
      soft_clip_back = exon_end - last_interval->end;
      return;

    // First/middle exon, exon end extends beyond guide interval
    } else if (last_interval->end < exon_end) {
      return;

    // Intervals don't support any TIDs
    } else if (exon_tids.empty()) {
      return;
    }

    // *^*^*^ *^*^*^ *^*^*^ *^*^*^
    // Next, update match TIDs and positions for read
    // *^*^*^ *^*^*^ *^*^*^ *^*^*^

    if (is_first_exon) {
      // Initialize matches with TIDs and positions
      for (const auto &tid : exon_tids) {
        uint match_pos = get_match_pos(first_interval, tid, read_strand, g2t,
                                       exon_start, used_backwards_overhang);
        matches.insert(std::make_tuple(tid, match_pos));
      }
      exon_tids.clear();

    } else {
      auto it_match = matches.begin();
      while (it_match != matches.end()) {
        if (exon_tids.count(std::get<0>(*it_match))) {
          ++it_match;
        } else {
          it_match = matches.erase(it_match);
        }
      }
      exon_tids.clear();

      if (matches.empty())
        return;
    }
  }

  // If we reach here, read is valid
  for (read_id_t &i : group) {
    CReadAln *group_read = bundle->reads[i];

    ReadInfo *this_read = new ReadInfo();
    this_read->valid_read = true;
    this_read->matches = matches;
    this_read->brec = group_read->brec;
    this_read->read_index = i;
    this_read->read_size = group_read->len;
    this_read->nh_i = matches.size();

    // Record mate information
    this_read->is_paired = (group_read->brec->flags() & BAM_FPAIRED);
    this_read->is_reverse = (group_read->brec->flags() & BAM_FREVERSE);
    read_info[i] = this_read;
  }
  return;
}

/**
 * Create bam_info entries
 *
 * @param final_transcripts tids to keep for this read
 * @param read_transcripts tids found for read
 * @param mate_transcripts tids found for mate
 * @param read_positions match positions by tid for read
 * @param mate_positions match positions by tid for mate
 * @param bam_info
 * @param read_index index of read in bundle->reads
 * @param mate_index index of mate in bundle->reads
 * @param read_size length of read
 * @param mate_size length of mate
 */
void add_mate_info(const std::set<tid_t> &final_transcripts,
                   const std::set<tid_t> &read_transcripts,
                   const std::set<tid_t> &mate_transcripts,
                   const std::map<tid_t, uint> &read_positions,
                   const std::map<tid_t, uint> &mate_positions,
                   std::map<read_id_t, ReadInfo *> &read_info, 
                   std::map<bam_id_t, BamInfo *> &bam_info, 
                   uint read_index, uint mate_index, uint mate_case) {

  auto this_read = read_info[read_index];
  auto mate_read = read_info[mate_index];

  if (mate_case == 1) {

    for (const tid_t &tid : final_transcripts) {

      // Add this read pair + transcript to bam_info
      bam_id_t n = bam_info.size();
      auto this_pair = new BamInfo();
      this_pair->same_transcript = true; // true if both mates map to same transcript
      this_pair->valid_pair = true;
      
      // Read 1 information
      this_pair->read_index = read_index;
      this_pair->tid = tid;
      this_pair->pos = read_positions.at(tid);
      this_pair->nh = this_read->nh_i;
      this_pair->read_size = this_read->read_size;
      this_pair->brec = this_read->brec;
      this_pair->is_reverse = this_read->is_reverse;

      // Read 2 information
      this_pair->mate_index = mate_index;
      this_pair->mate_tid = tid;
      this_pair->mate_pos = mate_positions.at(tid);
      this_pair->mate_nh = mate_read->nh_i;
      this_pair->mate_size = mate_read->read_size;
      this_pair->mate_brec = mate_read->brec;
      this_pair->mate_is_reverse = mate_read->is_reverse;

      bam_info[n] = this_pair;
    }

  } else if (mate_case == 2) {

    // Add this read pair + transcript to bam_info
    bam_id_t n = bam_info.size();
    auto this_pair = new BamInfo();
    this_pair->valid_pair = true;

    for (const tid_t &tid : read_transcripts) {
      // Read 1 information
      this_pair->read_index = read_index;
      this_pair->tid = tid;
      this_pair->pos = read_positions.at(tid);
      this_pair->nh = this_read->nh_i;
      this_pair->read_size = this_read->read_size;
      this_pair->brec = this_read->brec;
      this_pair->is_reverse = this_read->is_reverse;
    }

    for (const tid_t &tid : mate_transcripts) {
      // Read 2 information
      this_pair->mate_index = mate_index;
      this_pair->mate_tid = tid;
      this_pair->mate_pos = mate_positions.at(tid);
      this_pair->mate_nh = mate_read->nh_i;
      this_pair->mate_size = mate_read->read_size;
      this_pair->mate_brec = mate_read->brec;
      this_pair->mate_is_reverse = mate_read->is_reverse;
    }

    bam_info[n] = this_pair;
  }
  // could add more cases here if they exist

}

/**
 * Update read matches based on final_transcripts
 *
 * @param read_info current read information
 * @param final_transcripts tids to keep for this read
 */
void update_read_matches(ReadInfo *read_info,
                         const std::set<tid_t> &final_transcripts) {

  std::set<std::tuple<tid_t, uint>> new_matches;

  for (const auto &match : read_info->matches) {
    const tid_t &tid = std::get<0>(match);
    if (final_transcripts.find(tid) != final_transcripts.end()) {
      new_matches.insert(match);
    }
  }

  read_info->matches = std::move(new_matches);
}

/**
 * Process mate relationships using pre-established pairing information
 *
 * @param bundle current bundle of reads
 * @param bam_info read_info for all reads in bundle
 * @param mate_map map between read_id and mate_info
 */
void process_mate_pairs(BundleData *bundle, 
                        std::map<read_id_t, ReadInfo *> &read_info,
                        std::map<bam_id_t, BamInfo *> &bam_info) {

  GList<CReadAln> &reads = bundle->reads;
  const int MAX_MAPPINGS_PER_READ = 200;

  // Pre-compute read transcript sets to avoid repeated extraction
  std::vector<std::set<read_id_t>> read_transcript_cache(reads.Count());
  std::vector<std::map<read_id_t, uint>> read_position_cache(reads.Count());
  std::vector<bool> read_valid(reads.Count(), false);

  // *^*^*^ *^*^*^ *^*^*^ *^*^*^
  // First pass: Cache transcript information for valid reads
  // *^*^*^ *^*^*^ *^*^*^ *^*^*^

  for (int i = 0; i < reads.Count(); i++) {

    auto read_info_it = read_info.find(i);
    if (read_info_it == read_info.end() || !read_info_it->second->valid_read ||
        !read_info_it->second->is_paired) {
      continue;
    }

    auto this_read = read_info_it->second;

    // Skip reads with too many mappings
    if (this_read->matches.size() > max_mappings_per_read) {
      continue;
    }

    read_valid[i] = true;

    // Cache transcript sets
    for (const auto &match : this_read->matches) {
      tid_t tid = std::get<0>(match);
      uint pos = std::get<1>(match);
      read_transcript_cache[i].insert(tid);
      read_position_cache[i][tid] = pos;
    }
  }

  // *^*^*^ *^*^*^ *^*^*^ *^*^*^
  // Second pass: Process mate pairs using cached data
  // *^*^*^ *^*^*^ *^*^*^ *^*^*^

  // Track processed pairs to avoid duplicates
  std::set<std::pair<read_id_t, read_id_t>> processed_pairs;

  for (int i = 0; i < reads.Count(); i++) {
    if (!read_valid[i])
      continue;

    CReadAln *read = reads[i];
    std::string read_name = read->brec->name();
    auto this_read = read_info[i];

    // Process all known mate relationships for this read
    for (int j = 0; j < read->pair_idx.Count(); j++) {
      int mate_index = read->pair_idx[j];

      // Validate mate index and check if mate is valid
      if (mate_index < 0 || mate_index >= reads.Count() ||
          !read_valid[mate_index]) {
        read_valid[i] = false;
        this_read->valid_read = false;
        continue;
      }

      // Avoid processing the same pair twice (both directions)
      std::pair<read_id_t, read_id_t> pair_key =
          std::make_pair(std::min(i, mate_index), std::max(i, mate_index));
      if (processed_pairs.find(pair_key) != processed_pairs.end()) {
        continue;
      }
      processed_pairs.insert(pair_key);

      CReadAln *mate = reads[mate_index];
      std::string mate_name = mate->brec->name();
      auto mate_read = read_info[mate_index];

      // Use cached transcript sets for fast intersection
      const auto &read_transcripts = read_transcript_cache[i];
      const auto &mate_transcripts = read_transcript_cache[mate_index];
      const auto &read_positions = read_position_cache[i];
      const auto &mate_positions = read_position_cache[mate_index];

      // Find common transcripts
      std::set<tid_t> common_transcripts;
      std::set_intersection(
          read_transcripts.begin(), read_transcripts.end(),
          mate_transcripts.begin(), mate_transcripts.end(),
          std::inserter(common_transcripts, common_transcripts.begin()));

      // Apply mate pairing logic with early exit
      std::set<tid_t> final_transcripts;

      // *^*^*^ *^*^*^ *^*^*^ *^*^*^
      // Mate pair cases
      // *^*^*^ *^*^*^ *^*^*^ *^*^*^

      uint mate_case;

      if (!common_transcripts.empty()) {
        // Case 1: Mates share some transcripts - keep only shared ones
        final_transcripts = std::move(common_transcripts);
        mate_case = 1;

      } else if (read_transcripts.size() == 1 && mate_transcripts.size() == 1) {
        // Case 2: Each mate maps to exactly one transcript, but different ones
        // - allow it
        final_transcripts.insert(*read_transcripts.begin());
        final_transcripts.insert(*mate_transcripts.begin());
        mate_case = 2;

      } else if (read_transcripts.size() == 1 && mate_transcripts.size() == 0) {
        // Case 3: Read 1 mapped to 1 transcript, Read 2 mapped to none
        continue;

      } else if (mate_transcripts.size() == 1 && read_transcripts.size() == 0) {
        // Case 4: Read 2 mapped to 1 transcript, Read 1 mapped to none
        continue;

      } else {
        // Case 5: No common transcripts and at least one mate maps to multiple
        // - skip this pair
        continue;
      }

      // ALLOW: if pair is not mapped, then don't discard (even for multiple transcripts)
      // DISALLOW: if pair is mapped but to outside bundle, discard
      // we will test both ways with simulated data to see if that should really be disallowed

      // Update both reads' matches
      update_read_matches(this_read, final_transcripts);
      update_read_matches(mate_read, final_transcripts);

      // Add mate information for each valid transcript
      add_mate_info(final_transcripts, read_transcripts, mate_transcripts,
                    read_positions, mate_positions, read_info, bam_info, 
                    i, mate_index, mate_case);
    }
  }
}

void free_read_data(AlnGroups* aln_groups, std::map<read_id_t, ReadInfo *> &read_info,
                    std::map<bam_id_t, BamInfo *> &bam_info) {
  delete aln_groups;
  for (const auto &pair : read_info) {
    auto info = std::get<1>(pair);
    delete info;
  }
  for (const auto &pair : bam_info) {
    auto info = std::get<1>(pair);
    delete info;
  }
}

/**
 * Converts reads from bundle to new transcriptome-coordinate reads
 *
 * @param bundle current bundle of reads
 * @param io BAM reader and writer
 */
void convert_reads(BundleData *bundle, BamIO *io) {
  // Use bundle guides to build g2t tree
  std::unique_ptr<g2tTree> g2t = make_g2t_tree(bundle, io);
  if (g2t == nullptr)
    return;

  // Create read groups
  AlnGroups *aln_groups = new AlnGroups();
  GList<CReadAln> &reads = bundle->reads;

  for (int n = 0; n < reads.Count(); n++) {
    aln_groups->Add(reads[n], n);
  }

  // Create read_info to store info for every read
  std::map<read_id_t, ReadInfo *> read_info;

  // Create bam_info to store info for every output line
  std::map<bam_id_t, BamInfo *> bam_info;

  // First pass: process reads
  auto groups = aln_groups->groups;
  for (uint i = 0; i < groups.size(); i++) {
    uint n = groups[i][0]; // index of first read in group

    // Process the first read from this group
    process_read_out(bundle, n, g2t.get(), read_info, groups[i]);
  }

  // Second pass: process mate pairs
  process_mate_pairs(bundle, read_info, bam_info);

  // Third pass: write to BAM
  write_to_bam(io, bam_info, g2t.get());

  // Free allocated structures
  free_read_data(aln_groups, read_info, bam_info);
}
