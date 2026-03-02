use bramble_rs::{Cigar, CigarOp};
use bramble_rs::update_cigar_for_test;
use noodles::sam::alignment::record::cigar::{op::Kind as CigarKind, Op as SamCigarOp};

fn build_ideal(ops: &[(u32, CigarOp)]) -> Cigar {
    let mut cigar = Cigar::default();
    for (len, op) in ops {
        cigar.add_operation(*len, *op);
    }
    cigar
}

fn cigar_kinds(cigar: &noodles::sam::alignment::record_buf::Cigar) -> Vec<(CigarKind, usize)> {
    cigar.as_ref().iter().map(|op| (op.kind(), op.len())).collect()
}

#[test]
fn update_cigar_override_softclip_to_match() {
    let real = vec![
        SamCigarOp::new(CigarKind::SoftClip, 2),
        SamCigarOp::new(CigarKind::Match, 8),
    ];
    let ideal = build_ideal(&[
        (2, CigarOp::MatchOverride),
        (8, CigarOp::Match),
    ]);

    let (cigar, nm) = update_cigar_for_test(real, ideal);
    assert_eq!(nm, 0);
    assert_eq!(cigar_kinds(&cigar), vec![(CigarKind::Match, 10)]);
}

#[test]
fn update_cigar_inserts_deletion_from_ideal() {
    let real = vec![SamCigarOp::new(CigarKind::Match, 10)];
    let ideal = build_ideal(&[
        (5, CigarOp::Match),
        (1, CigarOp::Del),
        (4, CigarOp::Match),
    ]);

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
    let ideal = build_ideal(&[
        (5, CigarOp::Match),
        (1, CigarOp::Ins),
        (5, CigarOp::Match),
    ]);

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

/// Adjacent I+D in the ideal cigar should be merged into M (they cancel out).
/// Ideal: 5M 3I 3D 5M  →  merge_indels collapses the I+D into 3M → effective 13M.
/// Real: 10M (shorter than the expanded ideal; extra M's come from ideal).
#[test]
fn update_cigar_adjacent_indel_in_ideal_merges_to_match() {
    let real = vec![SamCigarOp::new(CigarKind::Match, 10)];
    let ideal = build_ideal(&[
        (5, CigarOp::Match),
        (3, CigarOp::Ins),
        (3, CigarOp::Del),
        (5, CigarOp::Match),
    ]);

    let (cigar, nm) = update_cigar_for_test(real, ideal);
    assert_eq!(nm, 0);
    assert_eq!(cigar_kinds(&cigar), vec![(CigarKind::Match, 13)]);
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
    let ideal = build_ideal(&[
        (2, CigarOp::ClipOverride),
        (8, CigarOp::Match),
    ]);

    let (cigar, nm) = update_cigar_for_test(real, ideal);
    assert_eq!(nm, 0);
    assert_eq!(
        cigar_kinds(&cigar),
        vec![(CigarKind::SoftClip, 2), (CigarKind::Match, 8)],
    );
}
