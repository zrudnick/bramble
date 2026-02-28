mod annotation;
mod alignment;
mod bam_input;
mod cli;
mod header;
mod g2t;
mod fasta;
mod evaluate;
mod pipeline;
mod sw;
mod types;

use anyhow::Result;
use clap::Parser;
use mimalloc::MiMalloc;
use tracing_subscriber::EnvFilter;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

fn main() -> Result<()> {
    let args = cli::Args::parse();
    
    // Initialize tracing subscriber
    let filter = EnvFilter::try_from_default_env()
        .unwrap_or_else(|_| {
            if args.quiet {
                EnvFilter::new("warn")
            } else {
                EnvFilter::new("info")
            }
        });
    tracing_subscriber::fmt()
        .with_env_filter(filter)
        .init();
    
    let transcripts = annotation::load_transcripts(&args.guide_gff)?;
    let fasta = if let Some(path) = &args.genome_fasta {
        Some(fasta::FastaDb::load(path)?)
    } else {
        None
    };
    let mut bam = bam_input::open_bam(&args.in_bam)?;
    let g2t = g2t::build_g2t(&transcripts, &bam.refname_to_id, fasta.as_ref())?;
    let out_header = header::build_header(&transcripts);
    let stats = pipeline::run(&args, &out_header, &mut bam, &g2t, fasta.as_ref())?;
    tracing::info!(
        total_reads = stats.total_reads,
        unmapped_reads = stats.unmapped_reads,
        read_groups = stats.read_groups,
        total_exons = stats.total_exons,
        "bramble-rs: processing complete"
    );
    Ok(())
}
