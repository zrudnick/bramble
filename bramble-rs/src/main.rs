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

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

fn main() -> Result<()> {
    let args = cli::Args::parse();
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
    if !args.quiet {
        eprintln!(
            "bramble-rs: header written, g2t built, reads scanned (total={}, unmapped={}, groups={}, exons={})",
            stats.total_reads, stats.unmapped_reads, stats.read_groups, stats.total_exons
        );
    }
    Ok(())
}
