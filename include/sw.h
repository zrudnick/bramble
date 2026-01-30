
#pragma once
#include "evaluate.h"

namespace bramble {

  struct Cigar;

  enum Anchor {
    START = 1,  // alignment must be near start of seq2
    END = 2     // alignment must be near end of seq2
  };

  enum CustomCigarOps {
    BAM_CMATCH_OVERRIDE = 10,
    BAM_CDEL_OVERRIDE = 11,
    BAM_CINS_OVERRIDE = 12,
    BAM_CLIP_OVERRIDE = 13
  };

  struct AlignmentResult {
    int score;
    Cigar cigar;
    int start_j, end_j;
    int start_i, end_i;
    int pos;          // pos (only relative to this guide exon!)
      // need to add pos_start from this guide exon too

    // Quality metrics
    int matches;
    int mismatches;
    int insertions;
    int deletions;
    double error_rate;  // (mismatches + insertions + deletions) / seq1_length
    double accuracy;    // matches / seq1_length
  };

  AlignmentResult 
  smith_waterman(const std::string &seq1, const std::string &seq2, 
                Anchor anchor);
}