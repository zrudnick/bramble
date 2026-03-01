// alignment.rs is used by the binary (BAM I/O path); suppress dead-code warnings
// when compiling as a library-only target.
#![allow(dead_code)]
use anyhow::{anyhow, Result};
use noodles::bam;
use noodles::sam::alignment::record::cigar::op::Kind;

#[derive(Debug, Clone, Copy)]
pub struct Segment {
    pub start: u32,
    pub end: u32,
}

/// Extract exon segments from a spliced alignment.
///
/// CIGAR `N` operations split exons (splice junctions).
/// Coordinates are returned as 1-based, half-open [start, end) to match gclib.
pub fn extract_exons(record: &bam::Record) -> Result<Vec<Segment>> {
    let Some(start_pos) = record.alignment_start() else {
        return Ok(Vec::new());
    };
    let start_pos = start_pos.map_err(|e| anyhow!(e))?;
    let mut ref_pos = start_pos.get() as u32;
    let mut exon_start = ref_pos;
    let mut exons: Vec<Segment> = Vec::new();

    for result in record.cigar().iter() {
        let op = result?;
        let len = op.len() as u32;
        match op.kind() {
            // Reference-consuming ops
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Deletion => {
                ref_pos = ref_pos.saturating_add(len);
            }
            Kind::Skip => {
                // End current exon before the intron
                if ref_pos > exon_start {
                    exons.push(Segment {
                        start: exon_start,
                        end: ref_pos,
                    });
                }
                ref_pos = ref_pos.saturating_add(len);
                exon_start = ref_pos;
            }
            // Non-reference-consuming ops: insertion, soft clip, hard clip, pad
            _ => {}
        }
    }

    if ref_pos > exon_start {
        exons.push(Segment {
            start: exon_start,
            end: ref_pos,
        });
    }

    Ok(exons)
}
