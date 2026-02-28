//! bramble-rs: project spliced genomic alignments into transcriptome space.
//!
//! # Library usage
//!
//! ```no_run
//! use bramble_rs::{GenomicAlignment, ProjectionConfig, project_group};
//! use bramble_rs::g2t::{G2TTree, build_g2t};
//! use bramble_rs::annotation::load_transcripts;
//! use std::collections::HashMap;
//!
//! // Build the genome-to-transcriptome index from a GTF/GFF.
//! // let transcripts = load_transcripts(path_to_gtf)?;
//! // let refname_to_id: HashMap<String, usize> = /* from BAM header or minimap2 index */;
//! // let index = build_g2t(&transcripts, &refname_to_id, None)?;
//! //
//! // let config = ProjectionConfig { long_reads: false, use_fasta: false };
//! //
//! // // Construct alignments from whatever source (minimap2-rs, noodles, etc.)
//! // let alns: Vec<GenomicAlignment> = /* … */;
//! // let projected = project_group(&alns, &index, &config);
//! ```

// Internal modules — not part of the public API.
pub(crate) mod alignment;
pub(crate) mod bam_input;
pub(crate) mod cli;
pub(crate) mod evaluate;
pub(crate) mod fasta;
pub(crate) mod header;
pub(crate) mod pipeline;
pub(crate) mod sw;
pub(crate) mod types;

// Public modules — stable API surface.
pub mod annotation;
pub mod g2t;

// Flat re-exports for the most commonly used public types.
pub use api::{GenomicAlignment, ProjectedAlignment, ProjectionConfig, project_group};
pub use evaluate::ReadEvaluationConfig;
pub use g2t::G2TTree;

// Re-exports needed by integration tests in tests/.
#[doc(hidden)]
pub use evaluate::{Cigar, CigarOp};
#[doc(hidden)]
pub use pipeline::update_cigar_for_test;

mod api;
