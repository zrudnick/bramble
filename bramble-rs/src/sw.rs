use crate::evaluate::{Cigar, CigarOp};

#[derive(Debug, Clone, Copy)]
pub enum Anchor {
    Start,
    End,
}

#[derive(Debug, Clone)]
pub struct AlignmentResult {
    pub score: i32,
    pub cigar: Cigar,
    pub start_j: i32,
    pub end_j: i32,
    pub start_i: i32,
    pub end_i: i32,
    pub pos: i32,
    #[allow(dead_code)]
    pub matches: i32,
    #[allow(dead_code)]
    pub mismatches: i32,
    #[allow(dead_code)]
    pub insertions: i32,
    #[allow(dead_code)]
    pub deletions: i32,
    #[allow(dead_code)]
    pub error_rate: f64,
    pub accuracy: f64,
}

pub fn smith_waterman(seq1: &[u8], seq2: &[u8], anchor: Anchor) -> AlignmentResult {
    let m = seq1.len();
    let n = seq2.len();

    if m == 0 || n == 0 {
        return AlignmentResult {
            score: -1,
            cigar: Cigar::default(),
            start_j: 0,
            end_j: 0,
            start_i: 0,
            end_i: 0,
            pos: 0,
            matches: 0,
            mismatches: 0,
            insertions: 0,
            deletions: 0,
            error_rate: 0.0,
            accuracy: 0.0,
        };
    }

    let match_score: i32 = 2;
    let mismatch_score: i32 = -1;
    let gap_score: i32 = -1;

    let cols = n + 1;
    // Use flat 1D arrays to avoid Vec<Vec<>> allocation overhead and enable
    // unchecked indexing.  Traceback is packed alongside the score to improve
    // cache locality (score in bits 2.., traceback in bits 0..1).
    let total = (m + 1) * cols;
    let mut score: Vec<i32> = vec![0i32; total];
    // Traceback: 0=stop, 1=diag, 2=up, 3=left
    let mut tb: Vec<u8> = vec![0u8; total];

    let mut max_score = 0i32;
    let mut max_i = 0usize;
    let mut max_j = 0usize;

    // SAFETY: All indices are in bounds because i ∈ [1, m] and j ∈ [1, n],
    // so i*cols+j ∈ [cols, (m+1)*cols-1] which is within the allocated `total`.
    // (i-1)*cols+(j-1) is always ≥ 0.  seq1[i-1] and seq2[j-1] are in bounds
    // because i-1 ∈ [0, m-1] and j-1 ∈ [0, n-1].
    unsafe {
        for i in 1..=m {
            let row = i * cols;
            let prev_row = (i - 1) * cols;
            let s1_byte = *seq1.get_unchecked(i - 1);

            for j in 1..=n {
                let s2_byte = *seq2.get_unchecked(j - 1);
                let match_bonus = if s1_byte == s2_byte { match_score } else { mismatch_score };

                let diag = *score.get_unchecked(prev_row + j - 1) + match_bonus;
                let up = *score.get_unchecked(prev_row + j) + gap_score;
                let left = *score.get_unchecked(row + j - 1) + gap_score;

                let cell = 0.max(diag.max(up.max(left)));
                *score.get_unchecked_mut(row + j) = cell;

                let tb_val = if cell == 0 {
                    0 // stop
                } else if cell == diag {
                    1 // diag
                } else if cell == up {
                    2 // up
                } else {
                    3 // left
                };
                *tb.get_unchecked_mut(row + j) = tb_val;

                if cell > max_score {
                    max_score = cell;
                    max_i = i;
                    max_j = j;
                }
            }
        }
    }

    let mut i = max_i;
    let mut j = max_j;
    let end_i = max_i as i32 - 1;
    let end_j = max_j as i32 - 1;

    let mut matches = 0i32;
    let mut mismatches = 0i32;
    let mut insertions = 0i32;
    let mut deletions = 0i32;
    let mut cigar = Cigar::default();

    while i > 0 && j > 0 && score[i * cols + j] > 0 {
        match tb[i * cols + j] {
            1 => {
                // diagonal
                cigar.add_operation(1, CigarOp::MatchOverride);
                if seq1[i - 1] == seq2[j - 1] {
                    matches += 1;
                } else {
                    mismatches += 1;
                }
                i -= 1;
                j -= 1;
            }
            3 => {
                // left
                cigar.add_operation(1, CigarOp::DelOverride);
                deletions += 1;
                j -= 1;
            }
            2 => {
                // up
                cigar.add_operation(1, CigarOp::InsOverride);
                insertions += 1;
                i -= 1;
            }
            _ => break,
        }
    }

    let start_i = i as i32;
    let start_j = j as i32;
    cigar.reverse();

    let pos = match anchor {
        Anchor::Start => n as i32 - end_j - 1,
        Anchor::End => start_j,
    };

    let mut score_out = max_score;
    match anchor {
        Anchor::Start => {
            let query_tail = start_i;
            let ref_tail = start_j;
            if query_tail > 0 && ref_tail > 0 {
                score_out = -1;
            } else if query_tail > ref_tail {
                let ins = (query_tail - ref_tail) as u32;
                score_out += gap_score * ins as i32;
                cigar.prepend_operation(ins, CigarOp::InsOverride);
                insertions += ins as i32;
            } else if ref_tail > query_tail {
                let del = (ref_tail - query_tail) as u32;
                score_out += gap_score * del as i32;
                cigar.prepend_operation(del, CigarOp::DelOverride);
                deletions += del as i32;
            }
        }
        Anchor::End => {
            let query_tail = m as i32 - 1 - end_i;
            let ref_tail = n as i32 - 1 - end_j;
            if query_tail > 0 && ref_tail > 0 {
                score_out = -1;
            } else if query_tail > ref_tail {
                let ins = (query_tail - ref_tail) as u32;
                score_out += gap_score * ins as i32;
                cigar.add_operation(ins, CigarOp::InsOverride);
                insertions += ins as i32;
            } else if ref_tail > query_tail {
                let del = (ref_tail - query_tail) as u32;
                score_out += gap_score * del as i32;
                cigar.add_operation(del, CigarOp::DelOverride);
                deletions += del as i32;
            }
        }
    }

    let total_errors = mismatches + insertions + deletions;
    let error_rate = if m == 0 { 0.0 } else { total_errors as f64 / m as f64 };
    let accuracy = if m == 0 { 0.0 } else { matches as f64 / m as f64 };

    AlignmentResult {
        score: score_out,
        cigar,
        start_j,
        end_j,
        start_i,
        end_i,
        pos,
        matches,
        mismatches,
        insertions,
        deletions,
        error_rate,
        accuracy,
    }
}
