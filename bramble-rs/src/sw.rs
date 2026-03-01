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

    let match_score = 2;
    let mismatch_score = -1;
    let gap_score = -1;

    let mut score = vec![vec![0i32; n + 1]; m + 1];
    let mut traceback = vec![vec![b'X'; n + 1]; m + 1];

    let mut max_score = 0i32;
    let mut max_i = 0usize;
    let mut max_j = 0usize;

    let s1 = seq1;
    let s2 = seq2;

    for i in 1..=m {
        for j in 1..=n {
            let diag = score[i - 1][j - 1]
                + if s1[i - 1] == s2[j - 1] {
                    match_score
                } else {
                    mismatch_score
                };
            let up = score[i - 1][j] + gap_score;
            let left = score[i][j - 1] + gap_score;

            let cell = 0.max(diag.max(up.max(left)));
            score[i][j] = cell;

            if cell == 0 {
                traceback[i][j] = b'X';
            } else if cell == diag {
                traceback[i][j] = b'D';
            } else if cell == up {
                traceback[i][j] = b'U';
            } else {
                traceback[i][j] = b'L';
            }

            if cell > max_score {
                max_score = cell;
                max_i = i;
                max_j = j;
            }
        }
    }

    let mut i = max_i;
    let mut j = max_j;
    let end_i = max_i as i32 - 1;
    let end_j = max_j as i32 - 1;

    let mut matches = 0;
    let mut mismatches = 0;
    let mut insertions = 0;
    let mut deletions = 0;
    let mut cigar = Cigar::default();

    while i > 0 && j > 0 && score[i][j] > 0 {
        match traceback[i][j] {
            b'D' => {
                cigar.add_operation(1, CigarOp::MatchOverride);
                if s1[i - 1] == s2[j - 1] {
                    matches += 1;
                } else {
                    mismatches += 1;
                }
                i -= 1;
                j -= 1;
            }
            b'L' => {
                cigar.add_operation(1, CigarOp::DelOverride);
                deletions += 1;
                j -= 1;
            }
            b'U' => {
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
