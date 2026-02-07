use clap::Parser;
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

    /// Reference annotation to guide conversion (GTF/GFF)
    #[arg(short = 'G', long = "guide", value_name = "GTF/GFF")]
    pub guide_gff: PathBuf,

    /// Output BAM path
    #[arg(short = 'o', long = "out", value_name = "BAM")]
    pub out_bam: PathBuf,

    /// Number of threads (CPUs) to use
    #[arg(short = 'p', long = "threads", default_value_t = 1)]
    pub threads: u8,

    /// Alignments are from long reads
    #[arg(long)]
    pub long: bool,

    /// Suppress progress bar and set logging level to WARN
    #[arg(short = 'q', long)]
    pub quiet: bool,

    /// Allow non-deterministic output order (read groups remain contiguous)
    #[arg(long)]
    pub unordered: bool,

    /// Unordered mode: flush output after this many records
    #[arg(long, default_value_t = 1024)]
    pub unordered_flush_records: usize,

    /// Genome sequence FASTA (optional)
    #[arg(short = 'S', long = "genome", value_name = "FASTA")]
    pub genome_fasta: Option<PathBuf>,
}
