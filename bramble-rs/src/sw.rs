use crate::evaluate::{Cigar, CigarOp};
use ksw2rs::{Aligner, Extz2Input, KSW_EZ_EXTZ_ONLY, KSW_EZ_APPROX_MAX, KSW_EZ_APPROX_DROP};

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

/// DNA5 scoring matrix (match=1, mismatch=-4) — matches C++ bramble ksw2 parameters.
/// Row-major 5x5 for alphabet {A=0, C=1, G=2, T=3, N=4}.
static DNA5_MAT: [i8; 25] = [
     1, -4, -4, -4, 0,  // A
    -4,  1, -4, -4, 0,  // C
    -4, -4,  1, -4, 0,  // G
    -4, -4, -4,  1, 0,  // T
     0,  0,  0,  0, 0,  // N
];

/// Convert ASCII base to 0-4 encoding for ksw2.
#[inline(always)]
fn encode_base(b: u8) -> u8 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 4, // N / unknown
    }
}

/// Encode an ASCII sequence into 0-4 encoding in-place into a provided buffer.
fn encode_seq(ascii: &[u8], buf: &mut Vec<u8>) {
    buf.clear();
    buf.reserve(ascii.len());
    for &b in ascii {
        buf.push(encode_base(b));
    }
}

pub fn smith_waterman(aligner: &mut Aligner, seq1: &[u8], seq2: &[u8], anchor: Anchor) -> AlignmentResult {
    let m = seq1.len();
    let n = seq2.len();

    if m == 0 || n == 0 {
        return AlignmentResult {
            score: -1,
            cigar: Cigar::default(),
            start_j: 0, end_j: 0,
            start_i: 0, end_i: 0,
            pos: 0,
            matches: 0, mismatches: 0, insertions: 0, deletions: 0,
            error_rate: 0.0, accuracy: 0.0,
        };
    }

    // For ksw2 extension alignment:
    // - query = seq1 (the soft-clipped read bases)
    // - target = seq2 (the reference exon sequence)
    //
    // For Anchor::End (left-clip rescue): we reverse both sequences so ksw2
    // extension runs right-to-left, matching the C++ align_reversed() approach.
    let mut qbuf = Vec::with_capacity(m);
    let mut tbuf = Vec::with_capacity(n);

    match anchor {
        Anchor::End => {
            // Reverse both sequences for left-extension
            let mut rev_q: Vec<u8> = seq1.to_vec();
            rev_q.reverse();
            let mut rev_t: Vec<u8> = seq2.to_vec();
            rev_t.reverse();
            encode_seq(&rev_q, &mut qbuf);
            encode_seq(&rev_t, &mut tbuf);
        }
        Anchor::Start => {
            encode_seq(seq1, &mut qbuf);
            encode_seq(seq2, &mut tbuf);
        }
    }

    let input = Extz2Input {
        query: &qbuf,
        target: &tbuf,
        m: 5,
        mat: &DNA5_MAT,
        q: 4,      // gap open penalty (matches C++ bramble)
        e: 1,      // gap extend penalty (matches C++ bramble)
        w: -1,     // unlimited band width
        zdrop: 40, // z-drop threshold (matches C++ bramble)
        end_bonus: 0,
        flag: KSW_EZ_EXTZ_ONLY | KSW_EZ_APPROX_MAX | KSW_EZ_APPROX_DROP,
    };

    let ez = aligner.align(&input);

    // Convert ksw2 packed CIGAR to our Cigar format and compute stats.
    let mut cigar = Cigar::default();
    let mut matches = 0i32;
    let mut mismatches = 0i32;
    let mut insertions = 0i32;
    let mut deletions = 0i32;
    let mut query_consumed = 0i32;
    let mut target_consumed = 0i32;

    let cigar_slice = &ez.cigar;

    // For Anchor::End (reversed alignment), iterate CIGAR in reverse
    // to restore original-orientation order (matching C++ approach).
    let cigar_iter: Box<dyn Iterator<Item = &u32>> = match anchor {
        Anchor::End => Box::new(cigar_slice.iter().rev()),
        Anchor::Start => Box::new(cigar_slice.iter()),
    };

    for &packed in cigar_iter {
        let len = (packed >> 4) as i32;
        let op = packed & 0xf;
        match op {
            0 => {
                // M (match/mismatch)
                cigar.add_operation(len as u32, CigarOp::MatchOverride);
                // Count actual matches/mismatches by comparing sequences
                let qi = query_consumed as usize;
                let ti = target_consumed as usize;
                for k in 0..len as usize {
                    let qpos = match anchor {
                        Anchor::End => m - 1 - (qi + k),
                        Anchor::Start => qi + k,
                    };
                    let tpos = match anchor {
                        Anchor::End => n - 1 - (ti + k),
                        Anchor::Start => ti + k,
                    };
                    if qpos < m && tpos < n && seq1[qpos] == seq2[tpos] {
                        matches += 1;
                    } else {
                        mismatches += 1;
                    }
                }
                query_consumed += len;
                target_consumed += len;
            }
            1 => {
                // I (insertion in query)
                cigar.add_operation(len as u32, CigarOp::InsOverride);
                insertions += len;
                query_consumed += len;
            }
            2 => {
                // D (deletion from query)
                cigar.add_operation(len as u32, CigarOp::DelOverride);
                deletions += len;
                target_consumed += len;
            }
            _ => {}
        }
    }

    // Compute coordinates in original (non-reversed) space.
    let (start_i, end_i, start_j, end_j, pos);
    match anchor {
        Anchor::Start => {
            // Extension from left: alignment starts at (0, 0) conceptually
            start_i = 0;
            end_i = query_consumed - 1;
            start_j = 0;
            end_j = target_consumed - 1;
            pos = start_j;
        }
        Anchor::End => {
            // Reversed extension: map back to original coordinates
            end_i = m as i32 - 1;
            start_i = end_i - query_consumed + 1;
            end_j = n as i32 - 1;
            start_j = end_j - target_consumed + 1;
            pos = n as i32 - end_j - 1;
        }
    }

    let score = ez.max as i32;
    let total_errors = mismatches + insertions + deletions;
    let error_rate = if m == 0 { 0.0 } else { total_errors as f64 / m as f64 };
    let accuracy = if m == 0 { 0.0 } else { matches as f64 / m as f64 };

    AlignmentResult {
        score,
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
