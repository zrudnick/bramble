pub type Tid = u32;
pub type ReadId = u32;
pub type RefId = i32;

// Fast hash maps / sets using AHash instead of the default SipHash.
// Import these throughout the codebase with `use crate::types::{HashMap, HashSet}`.
// Also import `HashMapExt` / `HashSetExt` when you need `::new()` or `::with_capacity()`.
pub(crate) type HashMap<K, V> = ahash::HashMap<K, V>;
pub(crate) type HashSet<K> = ahash::HashSet<K>;
pub(crate) use ahash::HashMapExt;
pub(crate) use ahash::HashSetExt;
