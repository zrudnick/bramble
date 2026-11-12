# bramble example data

A small data set for trying out bramble: simulated reads from transcripts on human
chromosome 22, with the annotation they were simulated from. The guides on the
[bramble website](https://zrudnick.github.io/bramble/) walk through running it:
[short reads](https://zrudnick.github.io/bramble/guide-short-reads.html) and
[long reads](https://zrudnick.github.io/bramble/guide-long-reads.html).

| File | Contents |
|---|---|
| `short_reads.fq.gz` | 99,777 single-end, forward-stranded 101-base reads (Illumina-like, simulated with Polyester) |
| `pacbio_reads.fa.gz` | 5,006 PacBio reads (simulated with PBSIM3 through SQANTI-SIM) |
| `ont_reads.fq.gz` | 5,163 Oxford Nanopore reads (simulated through SQANTI-SIM) |
| `chess2.2_chr22.gtf.gz` | chr22 of the expanded CHESS 2.2 annotation the reads were simulated from |

The genome sequence is not included; download chr22 of hg38 from UCSC:

```
curl -LO https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz
gunzip chr22.fa.gz
```

The reads are random subsets of the chr22 reads in one sample of each type from
the bramble paper's simulated data (short-read sample t1_s6, PacBio t2_s4,
ONT t0_s5), available at ftp.ccb.jhu.edu/pub/data/bramble.
