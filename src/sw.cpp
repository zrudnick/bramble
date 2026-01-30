
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

#include "bramble.h"
#include "sw.h"
#include "evaluate.h"

namespace bramble {

  AlignmentResult smith_waterman(const std::string &seq1, const std::string &seq2, 
                                Anchor anchor) {
    int m = seq1.length();
    int n = seq2.length();
    
    int match = 2;
    int mismatch = -1;
    int gap = -1;

    std::vector<std::vector<int>> score(m + 1, std::vector<int>(n + 1, 0));
    std::vector<std::vector<char>> traceback(m + 1, std::vector<char>(n + 1, 'X'));
    
    int max_score = 0;
    int max_i = 0, max_j = 0;
    
    for (int i = 1; i <= m; i++) {
      for (int j = 1; j <= n; j++) {
        
        int diag = score[i-1][j-1] + (seq1[i-1] == seq2[j-1] ? match : mismatch);
        int up = score[i-1][j] + gap;
        int left = score[i][j-1] + gap;
        
        score[i][j] = std::max({0, diag, up, left}); // include 0 for local alignment
        
        if (score[i][j] == 0) {
          traceback[i][j] = 'X';
        } else if (score[i][j] == diag) {
          traceback[i][j] = 'D';  // diagonal
        } else if (score[i][j] == up) {
          traceback[i][j] = 'U';  // up (deletion from seq2)
        } else {
          traceback[i][j] = 'L';  // left (insertion to seq2)
        }
        
        if (score[i][j] > max_score) {
          max_score = score[i][j];
          max_i = i;
          max_j = j;
        }
      }
    }
    
    int i = max_i, j = max_j;
    int end_i = max_i - 1, end_j = max_j - 1;
    
    int matches = 0;
    int mismatches = 0;
    int insertions = 0;
    int deletions = 0;
    
    Cigar cigar;
    while (i > 0 && j > 0 && score[i][j] > 0) {
      if (traceback[i][j] == 'D') {
        cigar.add_operation(1, BAM_CMATCH_OVERRIDE);
        if (seq1[i-1] == seq2[j-1]) matches++;
        else mismatches++;
        i--; j--;
      } else if (traceback[i][j] == 'L') {
        cigar.add_operation(1, BAM_CDEL_OVERRIDE);
        deletions++;
        j--;
      } else if (traceback[i][j] == 'U') {
        cigar.add_operation(1, BAM_CINS_OVERRIDE);
        insertions++;
        i--;
      } else break;
    }

    int start_i = i, start_j = j;
    cigar.reverse(); // traceback happens backwards

    int pos = 0;
    if (anchor == START) {
      pos = n - end_j - 1;
    } else if (anchor == END) {
      pos = start_j;
    }
  
    // Add operations to ends of local alignment
    if (anchor == START) {
      int query_tail = start_i;
      int ref_tail = start_j;

      if (query_tail > 0 && ref_tail > 0) {
        max_score = -1; // can't have I and D next to each other
      } 
      else if (query_tail > ref_tail) {
        int ins = query_tail - ref_tail;
        max_score += gap * ins;
        cigar.prepend_operation(ins, BAM_CINS_OVERRIDE);
        insertions += ins;
      }
      else if (ref_tail > query_tail) {
        int del = ref_tail - query_tail;
        max_score += gap * del;
        cigar.prepend_operation(del, BAM_CDEL_OVERRIDE);
        deletions += del;
      }
      
    } else if (anchor == END) {
      int query_tail = m - 1 - end_i;
      int ref_tail = n - 1 - end_j;

      if (query_tail > 0 && ref_tail > 0) {
        max_score = -1; // can't have I and D next to each other
      } 
      else if (query_tail > ref_tail) {
        int ins = query_tail - ref_tail;
        max_score += gap * ins;
        cigar.add_operation(ins, BAM_CINS_OVERRIDE);
        insertions += ins;
      }
      else if (ref_tail > query_tail) {
        int del = ref_tail - query_tail;
        max_score += gap * del;
        cigar.add_operation(del, BAM_CDEL_OVERRIDE);
        deletions += del;
      }
    }
    
    int total_errors = mismatches + insertions + deletions;
    double error_rate = (double)total_errors / m;
    double accuracy = (double)matches / m;
    
    return {max_score, cigar, start_j, end_j, start_i, end_i, pos, 
            matches, mismatches, insertions, deletions, error_rate, 
            accuracy};
  }
}