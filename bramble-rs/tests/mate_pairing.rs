//! End-to-end: when an aligner reports several alignments of one read pair at
//! the same positions (STAR does), the mates are matched by their HI tag, so
//! each alignment projects as a proper pair onto the transcript it fits.
//!
//! Annotation: tx1 is one exon chr1:1001-1400; tx2 has exons 1001-1100 and
//! 1301-1400. HI:1 is an unspliced pair (150M, fits tx1 only), HI:2 a spliced
//! pair (100M200N50M, fits tx2 only); all four records start at 1001.

use bramble_rs::{
    GenomicAlignment, ProjectionConfig, ProjectionContext, annotation::load_transcripts,
    g2t::build_g2t_from_refnames, project_group_with,
};
use std::io::Write;

fn mate(read1: bool, hi: i32, cigar: Vec<(u32, u8)>) -> GenomicAlignment {
    GenomicAlignment {
        query_name: "frag".to_string(),
        ref_id: 0,
        ref_start: 1001,
        is_reverse: !read1,
        cigar,
        sequence: None,
        is_paired: true,
        is_first_in_pair: read1,
        xs_strand: None,
        ts_strand: None,
        hit_index: hi,
        mate_ref_id: Some(0),
        mate_ref_start: Some(1001),
        mate_is_unmapped: false,
        read_len: 150,
    }
}

fn annotation() -> bramble_rs::G2TTree {
    let dir = std::env::temp_dir().join("bramble_mate_pairing_test");
    std::fs::create_dir_all(&dir).unwrap();
    let gtf = dir.join("anno.gtf");
    let mut f = std::fs::File::create(&gtf).unwrap();
    let rows = [
        ("transcript", 1001, 1400, "tx1"),
        ("exon", 1001, 1400, "tx1"),
        ("transcript", 1001, 1400, "tx2"),
        ("exon", 1001, 1100, "tx2"),
        ("exon", 1301, 1400, "tx2"),
    ];
    for (ty, s, e, id) in rows {
        writeln!(
            f,
            "chr1\ttest\t{ty}\t{s}\t{e}\t.\t+\t.\ttranscript_id \"{id}\"; gene_id \"g\";"
        )
        .unwrap();
    }
    drop(f);
    let txs = load_transcripts(&gtf).unwrap();
    build_g2t_from_refnames(&txs, &["chr1".to_string()], None).unwrap()
}

/// (transcript name, read1 input index, read2 input index) of each proper pair.
fn proper_pairs(gas: &[GenomicAlignment]) -> Vec<(String, usize, usize)> {
    let g2t = annotation();
    let proj = project_group_with(
        gas,
        &g2t,
        &ProjectionConfig::short_read(),
        &mut ProjectionContext::new(),
    );
    let mut pairs = Vec::new();
    for r1 in proj.iter().filter(|p| gas[p.input_index].is_first_in_pair) {
        for r2 in proj.iter().filter(|p| !gas[p.input_index].is_first_in_pair) {
            if r1.transcript_id == r2.transcript_id
                && r1.same_transcript_as_mate
                && r2.same_transcript_as_mate
                && r1.insert_size != 0
            {
                let name = g2t.transcript_name(r1.transcript_id).unwrap().to_string();
                pairs.push((name, r1.input_index, r2.input_index));
            }
        }
    }
    pairs.sort();
    pairs
}

#[test]
fn alignments_at_one_position_pair_by_hi() {
    let unspliced = || vec![(150, 0)];
    let spliced = || vec![(100, 0), (200, 3), (50, 0)];
    // the order STAR writes them in: read1 HI:2, read1 HI:1, read2 HI:1, read2 HI:2
    let gas = vec![
        mate(true, 2, spliced()),
        mate(true, 1, unspliced()),
        mate(false, 1, unspliced()),
        mate(false, 2, spliced()),
    ];
    // position-only pairing would join read1 HI:2 with read2 HI:1 (spliced with
    // unspliced: no shared transcript) and lose both pairs
    assert_eq!(
        proper_pairs(&gas),
        vec![("tx1".to_string(), 1, 2), ("tx2".to_string(), 0, 3)]
    );
}

#[test]
fn without_hi_a_lone_pair_still_pairs_by_position() {
    let gas = vec![
        mate(true, 0, vec![(150, 0)]),
        mate(false, 0, vec![(150, 0)]),
    ];
    assert_eq!(proper_pairs(&gas), vec![("tx1".to_string(), 0, 1)]);
}
