# Plan: Update and Reorganize bramble-rs

## Context

The `bramble-rs` Rust implementation was derived from the C++ bramble codebase but has drifted. The goals are:
1. **Output parity** with C++ — the projected alignments must match exactly.
2. **Dual-use**: work as both a standalone binary and a Rust library, callable from tools like `oarfish` and `minimap2-rs` without requiring BAM I/O.
3. **Idiomatic, performant Rust** — not a line-for-line port, but competitive with C++ throughput.
4. **Clean, testable codebase** with well-documented public API.

### Key findings from audit

- All 3 mate-pairing cases (unpaired / same-transcript / different-transcript) are already implemented.
- Default short-read config diverges from C++: Rust uses `max_clip=25, max_ins=5`; C++ uses `max_clip=2, max_ins=0` (stricter). Must fix.
- `bio = "3.0"` dependency is unused — remove.
- `lib.rs` exposes every module as `pub mod` with no curated public API — fix.
- NH computation for genomic multimappers was updated in C++ (`d925ddf`) and the Rust code needs to be synced.
- The clip-rescue Smith-Waterman is plain O(m×n) vs KSW2 in C++ — see note below.
- minimap2-rs `Mapping` and noodles `RecordBuf` are the two upstream input types to target.

---

## Phased Implementation Plan

### Phase 1 — Project structure: library + binary split

**Files to change:**
- `src/lib.rs` — replace current "pub everything" with a curated, documented public API surface
- `src/main.rs` — becomes a thin wrapper that constructs `G2TIndex`, reads BAM, calls library, writes BAM

**New file:** `src/api.rs`

Contains the library-facing types:

```rust
/// Strand of the alignment relative to the reference.
pub enum Strand { Forward, Reverse }

/// A genomic alignment record ready for projection.
/// Consumers construct this from noodles records or minimap2-rs Mappings.
pub struct GenomicAlignment {
    pub query_name: Arc<str>,
    pub ref_id: i32,
    pub ref_start: i64,
    pub strand: Strand,
    pub cigar: Vec<CigarOp>,
    pub sequence: Option<Vec<u8>>,
    pub is_paired: bool,
    pub is_first_in_pair: bool,
    pub xs_strand: Option<Strand>,   // XS tag (library strand for short reads)
    pub ts_strand: Option<Strand>,   // ts tag (long reads)
    pub hit_index: i32,              // HI tag value from input
    pub is_long_read: bool,
}

/// Result of projecting one alignment into transcriptome space.
pub struct ProjectedAlignment {
    pub transcript_id: u32,
    pub transcript_start: i64,
    pub strand: Strand,
    pub cigar: Vec<CigarOp>,
    pub similarity_score: f64,
    pub nh: u32,               // total hits for this read
    pub hi: u32,               // 1-based index of this hit
    pub is_primary: bool,
    pub same_transcript_as_mate: bool,
    pub insert_size: i64,
}
```

Provide `From<&noodles_sam::alignment::RecordBuf>` for `GenomicAlignment` so the binary can
convert with zero friction. Provide conversion hints (a helper module or doc example) for
minimap2-rs `Mapping`.

**Public API function (exposed from `lib.rs`):**

```rust
/// Project a group of alignments for the same query name.
/// Returns projected alignments sorted by (tid, hi).
pub fn project_group(
    alignments: &[GenomicAlignment],
    index: &G2TIndex,
    config: &ProjectionConfig,
) -> Vec<ProjectedAlignment>
```

Everything else (evaluate internals, pipeline threading, BAM writer) stays crate-private or
`pub(crate)`.

---

### Phase 2 — Fix configuration values and expose via CLI

**Files:** `src/evaluate.rs`, `src/cli.rs`

Update `ReadEvaluationConfig` defaults to match C++:

| Parameter | C++ short-read | C++ long-read | Current Rust | Fix to |
|-----------|---------------|--------------|-------------|--------|
| `max_clip` | 2 | 40 | 25 | Match C++ per mode |
| `max_ins` | 0 | 40 | 5 | Match C++ per mode |
| `max_gap` | — | 40 | 5 | Match C++ per mode |
| `similarity_threshold` | 0.80 | 0.60 | 0.80 | 0.60 for long reads |

Add a named constructor pair:
```rust
impl ReadEvaluationConfig {
    pub fn short_read() -> Self { /* C++ ShortReadEvaluator params */ }
    pub fn long_read() -> Self  { /* C++ LongReadEvaluator params */ }
}
```

Add CLI overrides in `src/cli.rs` (all optional, default to mode values):
- `--max-clip N`
- `--max-ins N`
- `--max-gap N`
- `--similarity F`

**`ProjectionConfig` (in `src/api.rs`):** wraps `ReadEvaluationConfig` and `use_fasta: bool`.

---

### Phase 3 — Fix NH/HI for genomic multimappers

**File:** `src/pipeline.rs`

The C++ commit `d925ddf` ("updated hit index to consider genomic multimappers") changed NH
computation to account for both read1 and read2 alignments when paired. Currently the Rust
code sets NH to `output_groups.len()` which only counts transcript matches for that read.

Fix: in `process_group_records()`, compute NH as the total count of projected alignments
across read AND mate (mirroring the C++ `total_matches` variable in `core.cpp` lines 238-256).

---

### Phase 4 — Dependency cleanup

**File:** `Cargo.toml`

- **Remove** `bio = "3.0"` — completely unused.
- **Replace** `crossfire = "3.0"` with `flume` (fast, battle-tested MPMC, widely used in the
  Rust ecosystem, semantics are simpler and better documented). `crossfire` is obscure and the
  `detect_backoff_cfg()` calls are non-standard.
- Keep `coitrees`, `noodles`, `ahash`, `mimalloc`, `tracing`, `clap`, `anyhow`, `indicatif`.

---

### Phase 5 — Library ergonomics and documentation

**Files:** `src/api.rs`, `src/lib.rs`

- Add `#[doc]` comments to all public types and the `project_group` function.
- Add a `pub use` re-export block in `lib.rs` so users only need `use bramble_rs::*` for the
  main types.
- Provide a usage example in `lib.rs` doc comment showing the minimap2-rs integration pattern:
  ```rust
  // let mappings: Vec<minimap2::Mapping> = aligner.map(...)?;
  // let records: Vec<GenomicAlignment> = mappings.iter().map(into_genomic_alignment).collect();
  // let projected = project_group(&records, &index, &config);
  ```

---

### Phase 6 — Testing

**Files:** `tests/`, new `tests/evaluate.rs`, `tests/integration.rs`

**Unit tests to add (`tests/evaluate.rs`):**
- `short_read_defaults_match_cpp` — construct `ReadEvaluationConfig::short_read()`, assert
  exact parameter values.
- `long_read_defaults_match_cpp` — same for long read.
- Expand `tests/update_cigar.rs` with:
  - Hard clip at boundary (preserved through merge)
  - Consecutive I+D interaction
  - Skip (`N`) operation stripped correctly
  - Override code interactions with soft clips (the ClipOverride → SoftClip path)

**Integration test (`tests/integration.rs`):**
- Test data will be provided locally on the development machine (too large to commit).
  A smaller fixture for CI can be extracted later.
- The test harness accepts the data path via env var (`BRAMBLE_TEST_DATA=/path/to/data`)
  and skips gracefully if the variable is not set, so `cargo test` still passes in CI.
- Run the Rust library projection and compare its BAM output against golden files generated
  from the C++ bramble binary.

---

## Execution Order

1. Phase 4 (dep cleanup — fast, no logic risk)
2. Phase 2 (config fix — fixes output parity immediately)
3. Phase 3 (NH/HI fix — output parity)
4. Phase 1 (library/binary split — architectural, builds on fixed logic)
5. Phase 5 (docs + ergonomics — polish)
6. Phase 6 (tests — validates everything)

---

## Critical files

| File | Role |
|------|------|
| `src/evaluate.rs` | Config defaults, segment building, similarity filter |
| `src/pipeline.rs` | NH/HI computation, threading, mate pairing |
| `src/g2t.rs` | Interval tree, overlap queries |
| `src/sw.rs` | Clip rescue alignment |
| `Cargo.toml` | Dependency changes |
| `src/lib.rs` | Public API surface |
| `src/api.rs` | New file: GenomicAlignment, ProjectedAlignment, ProjectionConfig |
| `src/types.rs` | Shared type aliases |
| `src/cli.rs` | New CLI flags |

---

## Out of scope (for now)

- Parallelizing the BAM writer beyond the current mutex approach.
- Supporting `--unordered` determinism guarantees beyond what the current implementation provides.

## Deferred but time-sensitive: clip rescue aligner

The current plain O(m×n) Smith-Waterman in `src/sw.rs` is functionally correct but may be a
meaningful throughput bottleneck for long-read workloads (clip rescue runs per read, on
sequences up to hundreds of bases). Rather than deferring this indefinitely, evaluate
performance impact early during Phase 6 integration testing. If the naive SW is a measured
bottleneck, replace it with one of:

1. **KSW2 FFI**: bind to the existing `subprojects/ksw2` C library already in the repo
   (already compiles on ARM64 via sse2neon). Add a `build.rs` to compile and link it,
   expose via a small `unsafe` wrapper. Highest performance, matches C++ exactly.
2. **`triple_accel`** crate: pure Rust SIMD-accelerated edit distance / alignment.
3. **`bio::alignment::pairwise`** banded DP: simpler, no FFI, reasonable performance.
