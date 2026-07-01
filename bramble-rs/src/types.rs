pub type Tid = u32;
pub type ReadId = u32;
pub type RefId = i32;

// Fast hash maps / sets using AHash instead of the default SipHash.
// Import these throughout the codebase with `use crate::types::{HashMap, HashSet}`.
// Also import `HashMapExt` / `HashSetExt` when you need `::new()` or `::with_capacity()`.
//
// DETERMINISM: we pin a FIXED hasher seed. AHash's own `RandomState::new()`
// (what `ahash::HashMapExt` uses) draws a per-process random seed, and `std`'s
// default `HashMap` is likewise randomly seeded — so iteration order varies
// run-to-run. bramble's projection output depends on map iteration order in two
// places: the annotation loader returns `transcripts.into_values()` (which fixes
// each transcript's dense `tid`), and the evaluator's per-read `matches`/`data`
// maps drive similarity tie-breaking. Random seeds therefore make the projected
// BAM/RAD (and any downstream quantification) non-reproducible. A fixed seed
// makes every map deterministic while keeping AHash's speed.

/// A deterministic fixed-seed AHash state. `with_seeds` is a const fn and does no
/// runtime RNG, so every map built from it hashes identically across runs.
#[doc(hidden)]
pub type DeterministicState = ahash::RandomState;

#[doc(hidden)]
#[inline]
pub fn deterministic_state() -> DeterministicState {
    // Arbitrary but fixed seeds; any constant works as long as it never changes.
    ahash::RandomState::with_seeds(
        0x243f_6a88_85a3_08d3,
        0x1319_8a2e_0370_7344,
        0xa409_3822_299f_31d0,
        0x082e_fa98_ec4e_6c89,
    )
}

#[doc(hidden)]
pub type HashMap<K, V> = std::collections::HashMap<K, V, DeterministicState>;
#[doc(hidden)]
pub type HashSet<K> = std::collections::HashSet<K, DeterministicState>;

/// Drop-in replacement for `ahash::HashMapExt`, but seeded deterministically.
#[doc(hidden)]
pub trait HashMapExt {
    fn new() -> Self;
    fn with_capacity(cap: usize) -> Self;
}
impl<K, V> HashMapExt for HashMap<K, V> {
    #[inline]
    fn new() -> Self {
        std::collections::HashMap::with_hasher(deterministic_state())
    }
    #[inline]
    fn with_capacity(cap: usize) -> Self {
        std::collections::HashMap::with_capacity_and_hasher(cap, deterministic_state())
    }
}

/// Drop-in replacement for `ahash::HashSetExt`, but seeded deterministically.
#[doc(hidden)]
pub trait HashSetExt {
    fn new() -> Self;
    fn with_capacity(cap: usize) -> Self;
}
impl<K> HashSetExt for HashSet<K> {
    #[inline]
    fn new() -> Self {
        std::collections::HashSet::with_hasher(deterministic_state())
    }
    #[inline]
    fn with_capacity(cap: usize) -> Self {
        std::collections::HashSet::with_capacity_and_hasher(cap, deterministic_state())
    }
}
