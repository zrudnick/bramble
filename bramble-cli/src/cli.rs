use clap::Parser;
use std::num::NonZeroU32;
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(
    name = "bramble-rs",
    about = "Project spliced genomic alignments into transcriptomic space",
    version
)]
pub struct Args {
    /// Input BAM with genomic alignments
    pub in_bam: PathBuf,

    /// Reference annotation for guiding alignment projection (GTF/GFF)
    #[arg(short = 'G', long = "guide", value_name = "GTF/GFF")]
    pub guide_gff: PathBuf,

    /// Output path/file name for the projected alignments
    #[arg(short = 'o', long = "out", value_name = "BAM")]
    pub out_bam: PathBuf,

    /// Number of threads (CPUs) to use
    #[arg(short = 'p', long = "threads", default_value_t = NonZeroU32::new(1).unwrap())]
    pub threads: NonZeroU32,

    /// Preset: alignments are from long reads
    #[arg(long)]
    pub lr: bool,

    /// Preset: alignments are from high-quality long reads
    #[arg(long = "lr-hq", conflicts_with = "lr")]
    pub lr_hq: bool,

    /// Deprecated spelling of --lr-hq; still accepted, but hidden from --help
    #[arg(long = "lr:hq", hide = true, conflicts_with = "lr")]
    pub lr_hq_deprecated: bool,

    /// Force strict boundary adherence
    #[arg(long)]
    pub strict: bool,

    /// Maximum added soft clip (overrides preset)
    #[arg(long = "max-soft-clip", value_name = "INT")]
    pub max_clip: Option<u32>,

    /// Maximum allowed insertion at any splice junction (overrides preset)
    #[arg(long = "max-junction-insertion", value_name = "INT")]
    pub max_junc_ins: Option<u32>,

    /// Maximum allowed deletion at any splice junction (overrides preset)
    #[arg(long = "max-junction-deletion", value_name = "INT")]
    pub max_junc_gap: Option<u32>,

    /// Candidate transcripts with similarity scores above this threshold are kept (overrides preset)
    #[arg(long, value_name = "FLOAT")]
    pub similarity_threshold: Option<f32>,

    /// Maximum allowed size of an insertion or deletion of an entire exon, for sequencing technologies that produce long indels (overrides preset)
    #[arg(long, value_name = "INT")]
    pub max_error_exon: Option<u32>,

    /// Suppress progress bar and set logging level to WARN
    #[arg(short = 'q', long)]
    pub quiet: bool,

    /// Allow non-deterministic output order (read groups remain contiguous)
    #[arg(long)]
    pub unordered: bool,

    /// Unordered mode: flush output after this many records
    #[arg(long, default_value_t = 1024)]
    pub unordered_flush_records: usize,

    /// Genome sequence file for long-read clip rescue (FASTA)
    #[arg(short = 'S', long = "genome", value_name = "FASTA")]
    pub genome_fasta: Option<PathBuf>,

    /// Assume stranded library, fr-secondstrand (read1 sense)
    #[arg(long, conflicts_with = "rf")]
    pub fr: bool,

    /// Assume stranded library, fr-firststrand, e.g. dUTP (read2 sense)
    #[arg(long, conflicts_with = "fr")]
    pub rf: bool,
}
