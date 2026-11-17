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
use rust_htslib::bam;
use rust_htslib::bam::Read as HtsRead;
use std::path::{Path, PathBuf};
use std::process::Command;

// ── helpers ──────────────────────────────────────────────────────────────────

fn test_data_dir() -> Option<PathBuf> {
    if let Ok(val) = std::env::var("BRAMBLE_TEST_DATA") {
        return Some(PathBuf::from(val));
    }
    let manifest = env!("CARGO_MANIFEST_DIR");
    let default_path = PathBuf::from(manifest).join("../test_data");
    if default_path.is_dir() {
        Some(default_path)
    } else {
        None
    }
}

fn count_bam_records(path: &Path) -> usize {
    let mut reader = bam::Reader::from_path(path).expect("open BAM");
    let mut count = 0usize;
    let mut record = bam::Record::new();
    loop {
        match reader.read(&mut record) {
            None => break,
            Some(Err(e)) => panic!("read_record error: {e}"),
            Some(Ok(())) => count += 1,
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
            eprintln!(
                "Skipping short_read_count_matches_golden: test data not found \
                       (set BRAMBLE_TEST_DATA or place test_data/ next to bramble-rs/)"
            );
            return;
        }
    };

    let gtf = data_dir.join("input/ref/chess2.2_ALL.gtf");
    let input = data_dir.join("input/short_reads/subsampled.bam");
    let golden = data_dir.join("output/short_reads/projected.bam");
    let out_path = std::env::temp_dir().join("bramble_rs_test_short.bam");

    let gtf_str = gtf.to_str().unwrap();
    let input_str = input.to_str().unwrap();
    run_binary(&[input_str, "-G", gtf_str, "-q"], &out_path);

    let got = count_bam_records(&out_path);
    let golden = count_bam_records(&golden);
    let _ = std::fs::remove_file(&out_path);

    assert_eq!(
        got, golden,
        "short-read record count mismatch: got {got}, expected {golden}"
    );
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

    let gtf = data_dir.join("input/ref/chess2.2_ALL.gtf");
    let input = data_dir.join("input/pacbio/long.0mm.gn.bam");
    let golden = data_dir.join("output/pacbio/long.0mm.tx.bam");
    let out_path = std::env::temp_dir().join("bramble_rs_test_long.bam");

    let gtf_str = gtf.to_str().unwrap();
    let input_str = input.to_str().unwrap();
    run_binary(&[input_str, "-G", gtf_str, "--long", "-q"], &out_path);

    let got = count_bam_records(&out_path);
    let golden = count_bam_records(&golden);
    let _ = std::fs::remove_file(&out_path);

    assert_eq!(
        got, golden,
        "long-read record count mismatch: got {got}, expected {golden}"
    );
}

// ── HI-aware mate pairing (self-contained fixture) ───────────────────────────

/// Write a 4-record SAM in the shape STAR reports one read pair with two
/// alignments at the same positions, told apart only by HI: HI:1 both mates
/// unspliced (150M), HI:2 both spliced (100M200N50M). Plus a GTF where the
/// unspliced pair fits tx1 (one exon 1001-1400) and the spliced pair fits tx2
/// (exons 1001-1100, 1301-1400).
fn write_hi_fixture(dir: &Path) -> (PathBuf, PathBuf) {
    use std::io::Write;
    std::fs::create_dir_all(dir).unwrap();
    let gtf = dir.join("anno.gtf");
    let mut f = std::fs::File::create(&gtf).unwrap();
    for (ty, s, e, id) in [
        ("transcript", 1001, 1400, "tx1"),
        ("exon", 1001, 1400, "tx1"),
        ("transcript", 1001, 1400, "tx2"),
        ("exon", 1001, 1100, "tx2"),
        ("exon", 1301, 1400, "tx2"),
    ] {
        writeln!(
            f,
            "chr1\ttest\t{ty}\t{s}\t{e}\t.\t+\t.\tgene_id \"g\"; transcript_id \"{id}\";"
        )
        .unwrap();
    }
    let sam = dir.join("pairs.sam");
    let mut f = std::fs::File::create(&sam).unwrap();
    writeln!(
        f,
        "@HD\tVN:1.6\tSO:unsorted\tGO:query\n@SQ\tSN:chr1\tLN:5000"
    )
    .unwrap();
    let seq = "ACGT".repeat(38)[..150].to_string();
    // STAR order: read1 HI:2, read1 HI:1, read2 HI:1, read2 HI:2
    for (flag, cigar, tlen, hi) in [
        (355, "100M200N50M", 350, 2),
        (99, "150M", 150, 1),
        (147, "150M", -150, 1),
        (403, "100M200N50M", -350, 2),
    ] {
        writeln!(
            f,
            "frag\t{flag}\tchr1\t1001\t3\t{cigar}\t=\t1001\t{tlen}\t{seq}\t*\tNH:i:2\tHI:i:{hi}\tAS:i:298"
        )
        .unwrap();
    }
    (sam, gtf)
}

/// (flag, transcript, pos, CIGAR, TLEN) of every output record, sorted.
fn bam_records(path: &Path) -> Vec<(u16, String, i64, String, i64)> {
    let mut reader = bam::Reader::from_path(path).expect("open BAM");
    let names: Vec<String> = reader
        .header()
        .target_names()
        .iter()
        .map(|n| String::from_utf8_lossy(n).into_owned())
        .collect();
    let mut out = Vec::new();
    for r in reader.records() {
        let r = r.expect("read record");
        out.push((
            r.flags(),
            names[r.tid() as usize].clone(),
            r.pos(),
            r.cigar().to_string(),
            r.insert_size(),
        ));
    }
    out.sort();
    out
}

/// Each alignment of the fragment projects as a proper pair onto the transcript
/// it fits. Pairing on positions alone would join read1 HI:2 with read2 HI:1
/// (spliced with unspliced: no shared transcript) and write nothing. If
/// `BRAMBLE_CPP` names the C++ binary, its output must be the same.
#[test]
fn alignments_at_one_position_pair_by_hi() {
    let dir = std::env::temp_dir().join("bramble_cli_hi_pairs");
    let (sam, gtf) = write_hi_fixture(&dir);
    let out = dir.join("rs.bam");
    run_binary(&[sam.to_str().unwrap(), "-G", gtf.to_str().unwrap()], &out);
    let got = bam_records(&out);
    let proper: Vec<_> = got.iter().filter(|r| r.0 & 0x2 != 0).collect();
    assert_eq!(proper.len(), 4, "records: {got:?}");
    let on = |tx: &str| proper.iter().filter(|r| r.1 == tx).count();
    assert_eq!((on("tx1"), on("tx2")), (2, 2), "records: {got:?}");

    if let Ok(cpp) = std::env::var("BRAMBLE_CPP") {
        let cpp_out = dir.join("cpp.bam");
        let status = Command::new(cpp)
            .args([sam.to_str().unwrap(), "-G", gtf.to_str().unwrap(), "-o"])
            .arg(&cpp_out)
            .status()
            .expect("failed to spawn C++ bramble");
        assert!(status.success(), "C++ bramble exited with status {status}");
        assert_eq!(bam_records(&cpp_out), got, "C++ and Rust outputs differ");
    }
}
