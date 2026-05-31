//! bramble-rs: project spliced genomic alignments into transcriptome space.
//!
//! # Library usage
//!
//! ```no_run
//! use bramble_rs::{GenomicAlignment, ProjectionConfig, project_group};
//! use bramble_rs::g2t::build_g2t_from_refnames;
//! use bramble_rs::annotation::load_transcripts;
//!
//! // Build the genome-to-transcriptome index from a GTF/GFF.
//! // let transcripts = load_transcripts(path_to_gtf)?;
//! // // `refnames[i]` is the chromosome whose 0-based RefId is `i`
//! // // (e.g. taken from a BAM header or a minimap2 target list).
//! // let refnames: Vec<String> = /* from BAM header or minimap2 index */;
//! // let index = build_g2t_from_refnames(&transcripts, &refnames, None)?;
//! //
//! // let config = ProjectionConfig { long_reads: false, use_fasta: false };
//! //
//! // // Construct alignments from whatever source (minimap2-rs, noodles, etc.)
//! // let alns: Vec<GenomicAlignment> = /* … */;
//! // let projected = project_group(&alns, &index, &config);
//! ```

// Internal modules. The BAM I/O + CLI (the old `bam_input`/`cli`/`header`/
// `pipeline` modules) now live in the separate `bramble-cli` binary crate so
// that this library does not depend on rust-htslib. The modules below are
// `#[doc(hidden)] pub` (not part of the stable API) only so `bramble-cli` can
// reach the projection internals it needs to build BAM records.
#[doc(hidden)]
pub mod alignment;
#[doc(hidden)]
pub mod cigar;
#[doc(hidden)]
pub mod evaluate;
#[doc(hidden)]
pub mod groups;
#[doc(hidden)]
pub mod types;
pub(crate) mod sw;

// Public modules — stable API surface.
pub mod annotation;
pub mod fasta;
pub mod g2t;

// Flat re-exports for the most commonly used public types.
pub use api::{
    GenomicAlignment, ProjectedAlignment, ProjectionConfig, ProjectionContext, project_group,
    project_group_with,
};
pub use evaluate::ReadEvaluationConfig;
pub use g2t::G2TTree;

// Re-exports used by bramble-cli and the integration tests.
#[doc(hidden)]
pub use evaluate::{Cigar, CigarOp};
#[doc(hidden)]
pub use cigar::update_cigar_for_test;

mod api;
