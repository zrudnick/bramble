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
