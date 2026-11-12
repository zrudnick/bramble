use bramble_rs::update_cigar_for_test;
use bramble_rs::{Cigar, CigarOp};
use noodles::sam::alignment::record::cigar::{Op as SamCigarOp, op::Kind as CigarKind};

fn build_ideal(ops: &[(u32, CigarOp)]) -> Cigar {
    let mut cigar = Cigar::default();
    for (len, op) in ops {
        cigar.add_operation(*len, *op);
    }
    cigar
}

fn cigar_kinds(cigar: &noodles::sam::alignment::record_buf::Cigar) -> Vec<(CigarKind, usize)> {
    cigar
        .as_ref()
        .iter()
        .map(|op| (op.kind(), op.len()))
        .collect()
}

#[test]
fn update_cigar_override_softclip_to_match() {
    let real = vec![
        SamCigarOp::new(CigarKind::SoftClip, 2),
        SamCigarOp::new(CigarKind::Match, 8),
    ];
    let ideal = build_ideal(&[(2, CigarOp::MatchOverride), (8, CigarOp::Match)]);

    let (cigar, nm) = update_cigar_for_test(real, ideal);
    assert_eq!(nm, 0);
    assert_eq!(cigar_kinds(&cigar), vec![(CigarKind::Match, 10)]);
}

#[test]
fn update_cigar_inserts_deletion_from_ideal() {
    let real = vec![SamCigarOp::new(CigarKind::Match, 10)];
    let ideal = build_ideal(&[(5, CigarOp::Match), (1, CigarOp::Del), (4, CigarOp::Match)]);

    let (cigar, nm) = update_cigar_for_test(real, ideal);
    assert_eq!(nm, 1);
    assert_eq!(
        cigar_kinds(&cigar),
        vec![
            (CigarKind::Match, 5),
            (CigarKind::Deletion, 1),
            (CigarKind::Match, 5),
        ]
    );
}

#[test]
fn update_cigar_inserts_insertion_from_ideal() {
    let real = vec![SamCigarOp::new(CigarKind::Match, 10)];
    let ideal = build_ideal(&[(5, CigarOp::Match), (1, CigarOp::Ins), (5, CigarOp::Match)]);

    let (cigar, nm) = update_cigar_for_test(real, ideal);
    assert_eq!(nm, 1);
    assert_eq!(
        cigar_kinds(&cigar),
        vec![
            (CigarKind::Match, 5),
            (CigarKind::Insertion, 1),
            (CigarKind::Match, 5),
        ]
    );
}

/// A leading hard-clip in the real CIGAR must survive the merge unchanged.
/// Real: 2H 10M  →  Ideal: 10M  →  Expected: 2H 10M
#[test]
fn update_cigar_preserves_leading_hard_clip() {
    let real = vec![
        SamCigarOp::new(CigarKind::HardClip, 2),
        SamCigarOp::new(CigarKind::Match, 10),
    ];
    let ideal = build_ideal(&[(10, CigarOp::Match)]);

    let (cigar, nm) = update_cigar_for_test(real, ideal);
    assert_eq!(nm, 0);
    assert_eq!(
        cigar_kinds(&cigar),
        vec![(CigarKind::HardClip, 2), (CigarKind::Match, 10)],
    );
}

/// Adjacent I+D in the ideal cigar are kept separate (matching C++ behavior).
/// Ideal: 5M 3I 3D 5M  →  output preserves I and D as-is.
/// Real: 10M (shorter than the expanded ideal; extra ops come from ideal).
#[test]
fn update_cigar_adjacent_indel_preserved() {
    let real = vec![SamCigarOp::new(CigarKind::Match, 10)];
    let ideal = build_ideal(&[
        (5, CigarOp::Match),
        (3, CigarOp::Ins),
        (3, CigarOp::Del),
        (5, CigarOp::Match),
    ]);

    let (cigar, nm) = update_cigar_for_test(real, ideal);
    assert_eq!(nm, 6);
    assert_eq!(
        cigar_kinds(&cigar),
        vec![
            (CigarKind::Match, 5),
            (CigarKind::Insertion, 3),
            (CigarKind::Deletion, 3),
            (CigarKind::Match, 5),
        ]
    );
}

/// A Skip (N / intron) operation in the real CIGAR is stripped during
/// `expand_real_cigar` and must not appear in the output.
/// Real: 5M 3N 5M  →  Ideal: 10M  →  Expected: 10M
#[test]
fn update_cigar_strips_skip_n_from_real() {
    let real = vec![
        SamCigarOp::new(CigarKind::Match, 5),
        SamCigarOp::new(CigarKind::Skip, 3),
        SamCigarOp::new(CigarKind::Match, 5),
    ];
    let ideal = build_ideal(&[(10, CigarOp::Match)]);

    let (cigar, nm) = update_cigar_for_test(real, ideal);
    assert_eq!(nm, 0);
    assert_eq!(cigar_kinds(&cigar), vec![(CigarKind::Match, 10)]);
}

/// `CigarOp::ClipOverride` (';') in the ideal forces the corresponding real
/// bases to SoftClip regardless of what the real CIGAR says.
/// Real: 10M  →  Ideal: 2 ClipOverride + 8 Match  →  Expected: 2S 8M
#[test]
fn update_cigar_clip_override_converts_match_to_softclip() {
    let real = vec![SamCigarOp::new(CigarKind::Match, 10)];
    let ideal = build_ideal(&[(2, CigarOp::ClipOverride), (8, CigarOp::Match)]);

    let (cigar, nm) = update_cigar_for_test(real, ideal);
    assert_eq!(nm, 0);
    assert_eq!(
        cigar_kinds(&cigar),
        vec![(CigarKind::SoftClip, 2), (CigarKind::Match, 8)],
    );
}

#[test]
fn update_cigar_insertion_in_clipped_overhang_joins_clip() {
    // Read 181M1I17M whose transcript ends after 172 bases: the insertion falls
    // inside the clipped overhang, so the whole overhang is one soft clip
    // (172M27S, as C++), not the invalid 172M9S1I17S.
    let real = vec![
        SamCigarOp::new(CigarKind::Match, 181),
        SamCigarOp::new(CigarKind::Insertion, 1),
        SamCigarOp::new(CigarKind::Match, 17),
    ];
    let ideal = build_ideal(&[(172, CigarOp::Match), (26, CigarOp::SoftClip)]);

    let (cigar, nm) = update_cigar_for_test(real, ideal);
    assert_eq!(nm, 0);
    assert_eq!(
        cigar_kinds(&cigar),
        vec![(CigarKind::Match, 172), (CigarKind::SoftClip, 27)]
    );
}

fn ops(v: &[(CigarKind, usize)]) -> Vec<SamCigarOp> {
    v.iter().map(|&(k, l)| SamCigarOp::new(k, l)).collect()
}

use CigarKind::{Deletion as D, Insertion as I, Match as M, SoftClip as S};
use bramble_rs::cigar::projected_nm;

#[test]
fn projected_nm_unchanged_alignment_keeps_mismatches() {
    let g = ops(&[(M, 150)]);
    assert_eq!(projected_nm(&g, &g, Some(1), Some(b"16C133"), 0, 0, 0), Some(1));
}

#[test]
fn projected_nm_insertion_absorbed_into_clip_is_not_an_edit() {
    // the CIGAR-bug read: 181M1I17M -> 172M27S
    let g = ops(&[(M, 181), (I, 1), (M, 17)]);
    let p = ops(&[(M, 172), (S, 27)]);
    assert_eq!(projected_nm(&g, &p, Some(1), Some(b"198"), 0, 0, 0), Some(0));
}

#[test]
fn projected_nm_mismatch_in_new_clip_dropped_with_md_kept_without() {
    let g = ops(&[(M, 100)]);
    let p = ops(&[(M, 90), (S, 10)]);
    // MD places the mismatch at position 95, inside the new clip
    assert_eq!(projected_nm(&g, &p, Some(1), Some(b"95A4"), 0, 0, 0), Some(0));
    // without MD the genome mismatch count is all we have
    assert_eq!(projected_nm(&g, &p, Some(1), None, 0, 0, 0), Some(1));
}

#[test]
fn projected_nm_counts_added_indels_and_md_deletions() {
    let g = ops(&[(M, 100)]);
    let p = ops(&[(M, 50), (D, 2), (M, 50)]);
    assert_eq!(projected_nm(&g, &p, Some(0), Some(b"100"), 0, 0, 0), Some(2));
    // genome deletion with a mismatch: MD 20A29^GG50 -> 1 mismatch + 2 deleted
    let g2 = ops(&[(M, 50), (D, 2), (M, 50)]);
    assert_eq!(projected_nm(&g2, &g2, Some(3), Some(b"20A29^GG50"), 0, 0, 0), Some(3));
}

#[test]
fn projected_nm_rescued_bases_use_rescue_mismatches() {
    // left 10 soft-clipped bases rescued as matches with 2 mismatches
    let g = ops(&[(S, 10), (M, 90)]);
    let p = ops(&[(M, 100)]);
    assert_eq!(projected_nm(&g, &p, Some(0), Some(b"90"), 10, 0, 2), Some(2));
}

#[test]
fn projected_nm_none_without_nm_or_md() {
    let g = ops(&[(M, 100)]);
    assert_eq!(projected_nm(&g, &g, None, None, 0, 0, 0), None);
}

#[test]
fn update_cigar_left_rescue_unrescued_bases_come_first() {
    // 22 soft-clipped bases, the leftward rescue reached the 16 next to the exon:
    // the 6 it did not reach stay clipped at the start (6S116M, as C++), not
    // between the rescue and the exon (invalid 16M6S100M)
    let real = vec![
        SamCigarOp::new(CigarKind::SoftClip, 22),
        SamCigarOp::new(CigarKind::Match, 100),
    ];
    let ideal = build_ideal(&[
        (6, CigarOp::ClipOverride),
        (16, CigarOp::MatchOverride),
        (100, CigarOp::Match),
    ]);
    let (cigar, _) = update_cigar_for_test(real, ideal);
    assert_eq!(
        cigar_kinds(&cigar),
        vec![(CigarKind::SoftClip, 6), (CigarKind::Match, 116)]
    );
}

#[test]
fn update_cigar_left_rescue_after_hard_clip() {
    // leading hard clip + soft clip: every soft-clipped base is matched against
    // the rescue ops (the padding loop used to stop front_hard bases early)
    let real = vec![
        SamCigarOp::new(CigarKind::HardClip, 5),
        SamCigarOp::new(CigarKind::SoftClip, 10),
        SamCigarOp::new(CigarKind::Match, 90),
    ];
    let ideal = build_ideal(&[
        (4, CigarOp::ClipOverride),
        (6, CigarOp::MatchOverride),
        (90, CigarOp::Match),
    ]);
    let (cigar, _) = update_cigar_for_test(real, ideal);
    assert_eq!(
        cigar_kinds(&cigar),
        vec![(CigarKind::HardClip, 5), (CigarKind::SoftClip, 4), (CigarKind::Match, 96)]
    );
}
