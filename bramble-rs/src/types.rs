pub type Tid = u32;
pub type ReadId = u32;
pub type RefId = i32;

// Fast hash maps / sets using AHash instead of the default SipHash.
// Import these throughout the codebase with `use crate::types::{HashMap, HashSet}`.
// Also import `HashMapExt` / `HashSetExt` when you need `::new()` or `::with_capacity()`.
#[doc(hidden)]
pub type HashMap<K, V> = ahash::HashMap<K, V>;
#[doc(hidden)]
pub type HashSet<K> = ahash::HashSet<K>;
#[doc(hidden)]
pub use ahash::HashMapExt;
#[doc(hidden)]
pub use ahash::HashSetExt;
