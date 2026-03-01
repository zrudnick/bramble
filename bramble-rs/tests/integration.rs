/// Integration tests comparing bramble-rs binary output against C++ golden BAMs.
///
/// These tests run the binary end-to-end and count projected records in the output.
/// They are skipped automatically when the test data is not present, so `cargo test`
/// still passes in CI without the test data checked out.
///
/// To run these tests locally:
///   cargo test --test integration
///
/// The test data is looked up from (in order):
///   1. The `BRAMBLE_TEST_DATA` environment variable.
///   2. The `../test_data` directory relative to this crate (the default location in the
///      mono-repo).
///
/// The golden BAMs were produced by the C++ bramble binary after fixing three long-read
/// bugs (correct_for_gaps UB, Rust get_strands_to_check, and C++ uint8_t refid overflow).
use noodles::bam;
use std::path::{Path, PathBuf};
use std::process::Command;

// ── helpers ──────────────────────────────────────────────────────────────────

fn test_data_dir() -> Option<PathBuf> {
    if let Ok(val) = std::env::var("BRAMBLE_TEST_DATA") {
        return Some(PathBuf::from(val));
    }
    let manifest = env!("CARGO_MANIFEST_DIR");
    let default_path = PathBuf::from(manifest).join("../test_data");
    if default_path.is_dir() { Some(default_path) } else { None }
}

fn count_bam_records(path: &Path) -> usize {
    let mut reader = bam::io::reader::Builder
        .build_from_path(path)
        .expect("open BAM");
    reader.read_header().expect("read header");
    let mut count = 0usize;
    let mut record = bam::Record::default();
    loop {
        match reader.read_record(&mut record) {
            Ok(0) => break,
            Ok(_) => count += 1,
            Err(e) => panic!("read_record error: {e}"),
        }
    }
    count
}

fn bramble_bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_bramble-rs"))
}

/// Run the binary, return the path to the output BAM.
fn run_binary(args: &[&str], out_path: &Path) {
    let status = Command::new(bramble_bin())
        .args(args)
        .arg("-o")
        .arg(out_path)
        .status()
        .expect("failed to spawn bramble-rs");
    assert!(status.success(), "bramble-rs exited with status {status}");
}

// ── tests ─────────────────────────────────────────────────────────────────────

/// Short-read projected record count must match the C++ golden BAM.
#[test]
fn short_read_count_matches_golden() {
    let data_dir = match test_data_dir() {
        Some(d) => d,
        None => {
            eprintln!("Skipping short_read_count_matches_golden: test data not found \
                       (set BRAMBLE_TEST_DATA or place test_data/ next to bramble-rs/)");
            return;
        }
    };

    let gtf      = data_dir.join("input/ref/chess2.2_ALL.gtf");
    let input    = data_dir.join("input/short_reads/subsampled.bam");
    let golden   = data_dir.join("output/short_reads/projected.bam");
    let out_path = std::env::temp_dir().join("bramble_rs_test_short.bam");

    let gtf_str   = gtf.to_str().unwrap();
    let input_str = input.to_str().unwrap();
    run_binary(&[input_str, "-G", gtf_str, "-q"], &out_path);

    let got    = count_bam_records(&out_path);
    let golden = count_bam_records(&golden);
    let _ = std::fs::remove_file(&out_path);

    assert_eq!(got, golden, "short-read record count mismatch: got {got}, expected {golden}");
}

/// Long-read projected record count must match the C++ golden BAM.
#[test]
fn long_read_count_matches_golden() {
    let data_dir = match test_data_dir() {
        Some(d) => d,
        None => {
            eprintln!("Skipping long_read_count_matches_golden: test data not found");
            return;
        }
    };

    let gtf      = data_dir.join("input/ref/chess2.2_ALL.gtf");
    let input    = data_dir.join("input/pacbio/long.0mm.gn.bam");
    let golden   = data_dir.join("output/pacbio/long.0mm.tx.bam");
    let out_path = std::env::temp_dir().join("bramble_rs_test_long.bam");

    let gtf_str   = gtf.to_str().unwrap();
    let input_str = input.to_str().unwrap();
    run_binary(&[input_str, "-G", gtf_str, "--long", "-q"], &out_path);

    let got    = count_bam_records(&out_path);
    let golden = count_bam_records(&golden);
    let _ = std::fs::remove_file(&out_path);

    assert_eq!(got, golden, "long-read record count mismatch: got {got}, expected {golden}");
}
