# bramble-cli

Command-line tool to project spliced **genomic** alignments into **transcriptomic** space.

`bramble-cli` is the command-line front-end for [Bramble](https://github.com/zrudnick/bramble).
It reads a genome-aligned BAM (with spliced, `N`-containing CIGARs), projects each alignment
onto compatible annotated transcripts, and writes a transcriptome-coordinate BAM ready for
quantifiers such as [Salmon](https://github.com/COMBINE-lab/salmon),
[oarfish](https://github.com/COMBINE-lab/oarfish), or TranSigner.

The projection logic lives in the I/O-free [`bramble-rs`](https://crates.io/crates/bramble-rs)
library; this crate adds the htslib-backed BAM reader/writer, threading, and CLI.

> **Note:** the installed executable is named **`bramble-rs`** (the crate is `bramble-cli` only
> to disambiguate it from the library on crates.io).

## Install

```sh
cargo install bramble-cli
```

This builds against the system [htslib](https://github.com/samtools/htslib) prerequisites
(`zlib`, `bzip2`); a C compiler and `libclang` (for `bindgen`) are required at build time.

## Usage

```sh
# Genome-aligned, name-collated BAM + annotation -> transcriptome-coordinate BAM.
# The input BAM is positional; -G/--guide is the annotation, -o/--out the output.
bramble-rs genome.bam -G annotation.gtf -o transcriptome.bam

# Long reads: enable long-read evaluation (--lr), and pass the genome FASTA
# (-S/--genome) to rescore soft-clipped ends.
bramble-rs genome.bam -G annotation.gtf -o transcriptome.bam \
    --lr -S genome.fa -p 16
```

Input BAMs should be **collated by read name** (e.g. `samtools collate`) so that all
alignments of a read are processed together. Run `bramble-rs --help` for the full set of
options (clip/gap tolerances, similarity threshold, strand handling, output ordering).

## License

MIT
