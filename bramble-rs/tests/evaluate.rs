use bramble_rs::ReadEvaluationConfig;

/// Verify that `ReadEvaluationConfig::short_read()` matches the C++ `ShortReadEvaluator`
/// constants (see `src/evaluate.cpp`, `ShortReadEvaluator::evaluate`).
#[test]
fn short_read_defaults_match_cpp() {
    let cfg = ReadEvaluationConfig::short_read();
    assert_eq!(cfg.max_clip, 5,       "max_clip");
    assert_eq!(cfg.max_ins, 0,        "max_ins");
    assert_eq!(cfg.max_gap, 0,        "max_gap");
    assert_eq!(cfg.max_junc_gap, 0,   "max_junc_gap");
    assert!((cfg.similarity_threshold - 0.90).abs() < 1e-6, "similarity_threshold");
    assert!(!cfg.ignore_small_exons,  "ignore_small_exons");
    assert_eq!(cfg.small_exon_size, 0, "small_exon_size");
}

/// Verify that `ReadEvaluationConfig::long_read()` matches the C++ `LongReadEvaluator`
/// constants (see `src/evaluate.cpp`, `LongReadEvaluator::evaluate`).
#[test]
fn long_read_defaults_match_cpp() {
    let cfg = ReadEvaluationConfig::long_read();
    assert_eq!(cfg.max_clip, 40,       "max_clip");
    assert_eq!(cfg.max_ins, 40,        "max_ins");
    assert_eq!(cfg.max_gap, 40,        "max_gap");
    assert_eq!(cfg.max_junc_gap, 40,   "max_junc_gap");
    assert!((cfg.similarity_threshold - 0.60).abs() < 1e-6, "similarity_threshold");
    assert!(cfg.ignore_small_exons,    "ignore_small_exons");
    assert_eq!(cfg.small_exon_size, 35, "small_exon_size");
}

/// `Default` should be equivalent to `short_read()`.
#[test]
fn default_is_short_read() {
    let default = ReadEvaluationConfig::default();
    let short   = ReadEvaluationConfig::short_read();
    assert_eq!(default.max_clip,            short.max_clip);
    assert_eq!(default.max_ins,             short.max_ins);
    assert_eq!(default.max_gap,             short.max_gap);
    assert_eq!(default.max_junc_gap,        short.max_junc_gap);
    assert_eq!(default.similarity_threshold, short.similarity_threshold);
    assert_eq!(default.ignore_small_exons,  short.ignore_small_exons);
    assert_eq!(default.small_exon_size,     short.small_exon_size);
}
