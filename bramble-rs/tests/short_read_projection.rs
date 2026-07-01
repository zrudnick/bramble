//! End-to-end regression for short-read projection: a single-exon (unspliced)
//! and a two-exon (spliced) short read must project onto their transcripts via
//! `project_group_with` with the default `ProjectionConfig::short_read()`.
//!
//! Guards the short-read `build_config` threshold bug: at similarity_threshold
//! 1.0 the `similarity > threshold` gate is unsatisfiable, so both reads
//! silently project 0 placements.

use bramble_rs::{
    GenomicAlignment, ProjectionConfig, ProjectionContext, annotation::load_transcripts,
    g2t::build_g2t_from_refnames, project_group_with,
};
use std::io::Write;

fn ga(name: &str, ref_start: i64, cigar: Vec<(u32, u8)>, read_len: usize) -> GenomicAlignment {
    GenomicAlignment {
        query_name: name.to_string(),
        ref_id: 0,
        ref_start,
        is_reverse: false,
        cigar,
        sequence: None,
        is_paired: false,
        is_first_in_pair: false,
        xs_strand: None,
        ts_strand: None,
        hit_index: 0,
        mate_ref_id: None,
        mate_ref_start: None,
        mate_is_unmapped: false,
        read_len,
    }
}

#[test]
fn salmon_fixture_projects() {
    let dir = std::env::temp_dir().join("bramble_salmon_fixture_diag");
    std::fs::create_dir_all(&dir).unwrap();
    let gtf = dir.join("anno.gtf");
    let mut f = std::fs::File::create(&gtf).unwrap();
    // tx1 single-exon chr1:101-300 (+). tx2 two-exon 500-600 + 800-900 (+).
    let rows = [
        ("chr1", "transcript", 101, 300, "+", "tx1"),
        ("chr1", "exon", 101, 300, "+", "tx1"),
        ("chr1", "transcript", 500, 900, "+", "tx2"),
        ("chr1", "exon", 500, 600, "+", "tx2"),
        ("chr1", "exon", 800, 900, "+", "tx2"),
    ];
    for (sn, ty, s, e, st, id) in rows {
        writeln!(
            f,
            "{sn}\ttest\t{ty}\t{s}\t{e}\t.\t{st}\t.\ttranscript_id \"{id}\"; gene_id \"g_{id}\";"
        )
        .unwrap();
    }
    drop(f);

    let txs = load_transcripts(&gtf).unwrap();
    eprintln!("loaded {} transcripts", txs.len());
    for t in &txs {
        eprintln!("  {} {} {} exons={:?}", t.id, t.seqname, t.strand, t.exons);
    }
    let refnames = vec!["chr1".to_string()];
    let g2t = build_g2t_from_refnames(&txs, &refnames, None).unwrap();
    eprintln!(
        "names={:?} lengths={:?}",
        g2t.transcript_names(),
        g2t.transcript_lengths()
    );

    let cfg = ProjectionConfig::short_read();
    let mut ctx = ProjectionContext::new();

    // Unspliced read fully inside tx1 exon.
    let r1 = ga("unspliced", 151, vec![(100, 0)], 100);
    let p1 = project_group_with(&[r1], &g2t, &cfg, &mut ctx);
    eprintln!("unspliced -> {} placements", p1.len());
    for p in &p1 {
        eprintln!(
            "  tid={} tstart={} tend={} qlen={} sim={} rev={}",
            p.transcript_id,
            p.transcript_start,
            p.transcript_end,
            p.query_aligned_len,
            p.similarity_score,
            p.is_reverse
        );
    }

    // Spliced read across tx2 junction: 551..600 (50M) | 199N | 800..849 (50M).
    let r2 = ga("spliced", 551, vec![(50, 0), (199, 3), (50, 0)], 100);
    let p2 = project_group_with(&[r2], &g2t, &cfg, &mut ctx);
    eprintln!("spliced -> {} placements", p2.len());
    for p in &p2 {
        eprintln!(
            "  tid={} tstart={} tend={} qlen={} sim={} rev={}",
            p.transcript_id,
            p.transcript_start,
            p.transcript_end,
            p.query_aligned_len,
            p.similarity_score,
            p.is_reverse
        );
    }

    assert!(!p1.is_empty(), "unspliced read must project onto tx1");
    assert!(!p2.is_empty(), "spliced read must project onto tx2");
}
