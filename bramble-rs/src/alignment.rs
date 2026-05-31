// Query exon segment type shared by the projection library and the bramble-cli
// binary. The htslib-based `extract_exons` (which reads a rust_htslib Record)
// lives in the bramble-cli crate, since the library does not depend on htslib.
#![allow(dead_code)]

/// A query exon segment: a maximal run of reference-consuming alignment between
/// splice junctions, as 1-based half-open `[start, end)` coordinates.
#[derive(Debug, Clone, Copy, Default)]
pub struct Segment {
    pub start: u32,
    pub end: u32,
}
