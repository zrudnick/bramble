# bramble🌿🫐 

bramble is a fast and versatile program that projects spliced genomic alignments into transcriptomic space, enabling transcript quantification to work with genomic alignments. It works with any alignment method and any quantification method that can use alignments as input. Furthermore, it has several modes optimized for different short- and long-read sequencing technologies. You can access bramble as a standalone executable in either C++ or Rust, as a library in Rust, or directly in Oarfish (and coming soon: Salmon!).

## Installation

### C++
For C++, bramble uses the `meson` build system. Currently, system dependencies are required for `zlib`, `bz2`, and `lzma`, though these could be built as dependencies later if wanted (an attempt will be made using wrapdb to build `zlib` and `bz2` if not found on the system). 

Also, to build you will need to have `meson`, `ninja`, and `autotools` installed.

To build you can either run build.sh directly, or do the following:

```
> meson setup build
```
then change to the build directory and build it

```
> cd build
> ninja
```

This should create the bramble executable!

This has been tested with meson 1.7 and 1.8, but likely also works with older versions.

### Rust
```
> cargo install bramble-cli
```

This builds against the system htslib prerequisites (zlib, bzip2); a C compiler and libclang (for bindgen) are required at build time.

## Quick start

### C++
```
# Genome-aligned, name-collated BAM + annotation -> transcriptome-coordinate BAM.
# The input BAM is positional; -G/--guide is the annotation, -o/--out the output.
bramble alignments.bam -G annotation.gtf -o transcriptome.bam

# Long reads: enable long-read evaluation (--lr), and pass the genome FASTA
# (-S/--genome) to recover soft-clipped ends.
bramble genome.bam -G annotation.gtf -o transcriptome.bam \
   --lr -S genome.fa -p 16
```

### Rust
```
# Genome-aligned, name-collated BAM + annotation -> transcriptome-coordinate BAM.
# The input BAM is positional; -G/--guide is the annotation, -o/--out the output.
bramble-rs genome.bam -G annotation.gtf -o transcriptome.bam

# Long reads: enable long-read evaluation (--lr), and pass the genome FASTA
# (-S/--genome) to recover soft-clipped ends.
bramble-rs genome.bam -G annotation.gtf -o transcriptome.bam \
    --lr -S genome.fa -p 16
```

## Usage

```
bramble <in.bam> [-G <annotation.gtf>] [-o <out.bam>] [-p <UINT>] [-S <genome.fa>]
 [--help] [--version] [--quiet] [--fr] [--rf] [--lr] [-lr-hq] [--strict]
 [--max-soft-clip <UINT>] [--max-junction-insertion <UINT>] [--max-junction-deletion <UINT>]
 [--max-error-exon <UINT>] [--similarity-threshold <FLOAT>]

Options:
 --help                   : print this usage message and exit
 --version                : print just the version at stdout and exit
 --quiet                  : turn off verbose mode (log bundle processing details)
 --fr                     : assume stranded library (first-strand, read2 sense)
 --rf                     : assume stranded library (second-strand, read1 sense)
 --lr                     : preset: alignments are from long reads
 --lr-hq                  : preset: alignments are from high-quality long reads
 --strict                 : force strict boundary adherence
 --max-soft-clip          : maximum added soft clip
 --max-junction-insertion : maximum allowed insertion at any splice junction
 --max-junction-deletion  : maximum allowed deletion at any splice junction
 --max-error-exon         : maximum allowed size of an insertion or deletion of an entire exon, for sequencing technologies that produce long indels
 --similarity-threshold   : candidate transcripts with similarity scores above this threshold are kept
 -G, --guide <file>       : reference annotation for guiding alignment projection (GTF/GFF)
 -S, --genome <file>      : genome sequence file for long-read clip rescue (FASTA format)\n\
 -o, --out <file>         : output path/file name for the projected alignments (default: stdout)
 -p <int>                 : number of threads (CPUs) to use (default: 1)
```

## Rust library

bramble-rs is the Rust library at the core of bramble. This crate is I/O-free: it has no dependency on htslib / rust-htslib. You bring the alignments (from rust-htslib, noodles, minimap2-rs, etc.) as plain structs, and bramble returns projected alignments as plain structs. The companion bramble-cli crate wraps this library in a BAM-in/BAM-out command-line tool (the bramble-rs executable).

See [bramble-rs](https://github.com/zrudnick/bramble/blob/main/bramble-rs/README.md) for usage examples.

## Getting Help
Run ```bramble --help``` to see all available options.
For issues or questions, please open an issue!

You may also contact us for questions or suggestions at zrudnic1@jh.edu.

## License
Bramble is free, open source software released under an [MIT License](https://opensource.org/license/MIT).
