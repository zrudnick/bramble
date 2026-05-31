// Binary-side alignment helpers: extract query exon segments from a rust_htslib
// BAM record. The `Segment` type itself lives in the bramble-rs library.
#![allow(dead_code)]
use anyhow::Result;
pub use bramble_rs::alignment::Segment;
use rust_htslib::bam::record::{Cigar, Record};

/// Extract exon segments from a spliced alignment.
///
/// CIGAR `N` operations split exons (splice junctions).
/// Coordinates are returned as 1-based, half-open [start, end) to match gclib.
pub fn extract_exons(record: &Record) -> Result<Vec<Segment>> {
    if record.is_unmapped() || record.pos() < 0 {
        return Ok(Vec::new());
    }
    // pos() is 0-based; convert to 1-based to match noodles convention used downstream.
    let mut ref_pos = (record.pos() + 1) as u32;
    let mut exon_start = ref_pos;
    let mut exons: Vec<Segment> = Vec::new();

    for op in record.cigar().iter() {
        match op {
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) | Cigar::Del(n) => {
                ref_pos = ref_pos.saturating_add(*n);
            }
            Cigar::RefSkip(n) => {
                if ref_pos > exon_start {
                    exons.push(Segment { start: exon_start, end: ref_pos });
                }
                ref_pos = ref_pos.saturating_add(*n);
                exon_start = ref_pos;
            }
            // Non-reference-consuming: Ins, SoftClip, HardClip, Pad
            _ => {}
        }
    }

    if ref_pos > exon_start {
        exons.push(Segment { start: exon_start, end: ref_pos });
    }

    Ok(exons)
}
