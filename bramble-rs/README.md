# bramble-rs

Project spliced **genomic** alignments into **transcriptomic** coordinates.

`bramble-rs` is the Rust library at the core of [Bramble](https://github.com/zrudnick/bramble).
It takes alignments of reads against a genome (with spliced, `N`-containing CIGARs) plus a
transcript annotation, and projects each alignment onto the transcripts it is compatible with —
re-deriving the transcriptome-space coordinates, CIGAR, and an exon-chain compatibility score.
This lets you align once against the genome (preserving the discovery power of spliced alignment)
and still feed transcriptome-coordinate alignments to quantifiers that expect them.

This crate is **I/O-free**: it has no dependency on htslib / `rust-htslib`. You bring the
alignments (from `rust-htslib`, [`noodles`](https://crates.io/crates/noodles),
[`minimap2-rs`](https://crates.io/crates/minimap2), etc.) as plain structs, and bramble returns
projected alignments as plain structs. The companion [`bramble-cli`](https://crates.io/crates/bramble-cli)
crate wraps this library in a BAM-in/BAM-out command-line tool (the `bramble-rs` executable).

## Example

```rust,no_run
use std::path::Path;
use bramble_rs::{GenomicAlignment, ProjectionConfig, project_group};
use bramble_rs::annotation::load_transcripts;
use bramble_rs::g2t::build_g2t_from_refnames;

# fn main() -> anyhow::Result<()> {
// 1. Build the genome→transcriptome index from a GTF/GFF.  `refnames[i]` is the
//    chromosome whose 0-based reference id is `i` (e.g. from a BAM header or a
//    minimap2 target list), so the index can map an alignment's `ref_id` to a chromosome.
let transcripts = load_transcripts(Path::new("annotation.gtf"))?;
let refnames: Vec<String> = vec![/* "chr1", "chr2", ... */];
let index = build_g2t_from_refnames(&transcripts, &refnames, None)?;

// 2. For each read, gather its genomic alignments (all alignments of one read are a
//    "group") and project them.  Construct `GenomicAlignment`s from any source.
let config = ProjectionConfig { long_reads: true, use_fasta: false, junc_miss_discount: 1.0 };
let alns: Vec<GenomicAlignment> = vec![/* ... */];
let projected = project_group(&alns, &index, &config);

for p in &projected {
    // p.transcript_id indexes into the transcripts used to build `index`.
    println!("{} -> tx {} @ {}..{} (score {})",
        "read", p.transcript_id, p.transcript_start, p.transcript_end, p.similarity_score);
}
# Ok(())
# }
```

For higher throughput, reuse a [`ProjectionContext`] across calls with
[`project_group_with`] to avoid re-allocating scratch buffers per read.

## Soft-clip rescue

When `ProjectionConfig::use_fasta` is set and the index is built with a genome FASTA,
bramble uses Smith–Waterman (via `ksw2`) to rescore soft-clipped ends against the candidate
transcript, recovering alignments that genome-aligner clipping would otherwise drop. This is
most useful for long reads.

## Key types

- [`GenomicAlignment`] — one input alignment of a read against the genome.
- [`ProjectionConfig`] — long-read mode, soft-clip-rescue toggle, junction-miss discount.
- [`ProjectedAlignment`] — one output alignment in transcript space (`transcript_id`,
  `transcript_start`/`transcript_end`, `aligned_len`, `query_aligned_len`, `similarity_score`, strand).
- [`G2TTree`] — the genome→transcriptome interval index; also exposes transcript names/lengths.
- [`project_group`] / [`project_group_with`] — project one read's alignment group.

## License

MIT
