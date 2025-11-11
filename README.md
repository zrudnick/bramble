# Bramble: project spliced genomic alignments into transcriptomic space🌿🫐 
Bramble converts genomic BAM alignments to transcriptomic coordinates, enabling accurate transcript quantification while preserving the discovery power of genome alignments.

## Why use Bramble?

Traditional RNA-seq quantification involves a fundamental trade-off: lightweight pseudoalignment and alignment-based methods may efficiently estimate transcript abundance, but they are limited by dependence on existing annotations. They miss novel isoforms and splice junctions, divergent sequences (SNPs, mutations, strain variants), unannotated transcripts, and transcriptional noise. Genome-based methods capture all of these features through spliced alignment, but quantification tools expect transcriptomic coordinates, creating a compatibility barrier.

Bramble bridges this gap by projecting spliced genomic alignments into transcriptomic space, giving you *the best of both worlds*: discover novel isoforms, variants, and unannotated features through genome alignment, and also run standard quantification tools that expect transcriptomic coordinates. Some benefits of using Bramble include:

*Improved quantification accuracy:* 
- More accurate downstream expression estimates compared to direct alignment-based transcriptome quantification
- Better handling of multimapped reads in transcript context
- Annotation-guided filtering removes noise and spurious alignments

*Flexible workflow:*
- Works with your existing genome aligner (HISAT2, STAR, minimap2, etc.)
- Compatible with downstream tools expecting transcriptomic BAMs (Salmon, Oarfish, TranSigner, etc.)
- Supports both short and long-read data

# Quick start

For the `meson` build system: currently, system dependencies are required for 
`zlib`, `bz2`, and `lzma`, though these could be built as dependencies later
if wanted (an attempt will be made using wrapdb to build `zlib` and `bz2` if
not found on the system).

Also, to build you will need to have `meson`, `ninja`, and `autotools` installed.

To build you can either run `build.sh` directly, or do the following:


```
> meson setup build
```

then change to the build directory and build it

```
> cd build
> ninja
```

This should create the `bramble` executable!

This has been tested with `meson` 1.7 and 1.8, but likely also works with older versions.

# Usage 📝 

*Basic command:*
```
bramble input.bam -G annotation.gtf -o output.bam
```

## Full options:

```
bramble <in.bam ...> [-G <guide_gff>] [-o <out.bam>] [-p <cpus>] [options]
```

*Required:*
```
<in.bam>: One or more input BAM files with genomic alignments
-G <file>: Reference annotation (GTF/GFF) to guide coordinate conversion
-o <file>: Output BAM file path
```
*Options:*
```
-p <int>: Number of threads to use (default: 1)
--long: Optimize for long-read alignments
--verbose: Print detailed processing information
--help: Show help message
--version: Print version and exit
```

# Getting help
Run ```bramble --help``` to see all available options.
For issues or questions, please open an issue!

You may also contact us for questions or suggestions at zrudnic1@jh.edu.

# License
