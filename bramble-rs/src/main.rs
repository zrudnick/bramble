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
use std::mem::ManuallyDrop;
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

    // Load annotation and FASTA in parallel — they are independent I/O-bound
    // tasks (~2.9s and ~3.4s respectively) so overlapping them saves ~2.9s.
    let (transcripts, fasta) = std::thread::scope(|s| {
        let gtf_handle = s.spawn(|| annotation::load_transcripts(&args.guide_gff));
        let fasta_handle = s.spawn(|| {
            args.genome_fasta
                .as_ref()
                .map(|path| fasta::FastaDb::load(path))
                .transpose()
        });
        let transcripts = gtf_handle.join().expect("annotation thread panicked")?;
        let fasta = fasta_handle.join().expect("fasta thread panicked")?;
        Ok::<_, anyhow::Error>((transcripts, fasta))
    })?;

    let mut bam = bam_input::open_bam(&args.in_bam)?;
    let g2t = g2t::build_g2t(&transcripts, &bam.refname_to_id, fasta.as_ref())?;
    let hts_header = header::build_hts_header(&transcripts);

    // Optimization: wrap large data structures in ManuallyDrop so their
    // destructors are skipped when main() returns.  The OS reclaims all
    // process memory on exit, so running destructors for the g2t interval
    // trees, FASTA database, and transcript list is pure overhead (~6% of
    // wall time in profiling).
    let g2t = ManuallyDrop::new(g2t);
    let fasta = ManuallyDrop::new(fasta);
    let _transcripts = ManuallyDrop::new(transcripts);

    let stats = pipeline::run(&args, &hts_header, &mut bam, &g2t, fasta.as_ref())?;
    tracing::info!(
        total_reads = stats.total_reads,
        unmapped_reads = stats.unmapped_reads,
        read_groups = stats.read_groups,
        total_exons = stats.total_exons,
        "bramble-rs: processing complete"
    );
    Ok(())
}
