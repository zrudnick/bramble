// header.rs is used only by the binary.
#![allow(dead_code)]
use bramble_rs::annotation::Transcript;
use rust_htslib::bam::Header;
use rust_htslib::bam::header::HeaderRecord;

// pub fn build_hts_header(transcripts: &[Transcript]) -> rust_htslib::bam::Header {
//     use rust_htslib::bam::header::HeaderRecord;
//     let mut header = rust_htslib::bam::Header::new();
//     let mut hd = HeaderRecord::new(b"HD");
//     hd.push_tag(b"VN", "1.0");
//     header.push_record(&hd);
//     for tx in transcripts {
//         let len = tx.length();
//         if len == 0 {
//             continue;
//         }
//         let mut sq = HeaderRecord::new(b"SQ");
//         sq.push_tag(b"SN", &tx.id);
//         sq.push_tag(b"LN", len);
//         header.push_record(&sq);
//     }
//     header
// }

pub fn build_hts_header(
    transcripts: &[Transcript],
    input_header: &Header,
    pg_args: &str,
) -> Header {
    let mut header = Header::new();

    // HD record
    // Records are grouped by read name, not sorted (and never by coordinate;
    // oarfish rejects SO:coordinate). Same line as the C++ version.
    let mut hd = HeaderRecord::new(b"HD");
    hd.push_tag(b"VN", "1.6");
    hd.push_tag(b"SO", "unsorted");
    hd.push_tag(b"GO", "query");
    header.push_record(&hd);

    // SQ records from transcripts
    for tx in transcripts {
        let len = tx.length();
        if len == 0 {
            continue;
        }
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", &tx.id);
        sq.push_tag(b"LN", len);
        header.push_record(&sq);
    }

    // Retain the PG and RG records from the input BAM header (RG: the read
    // group tags passed through on records must stay declared), as C++ does
    let header_view = input_header.to_hashmap();
    for kind in ["PG", "RG"] {
        if let Some(records) = header_view.get(kind) {
            for r in records {
                let mut rec = HeaderRecord::new(kind.as_bytes());
                for (key, val) in r {
                    rec.push_tag(key.as_bytes(), val);
                }
                header.push_record(&rec);
            }
        }
    }
    // CO lines have no tags; copy them verbatim
    for co in input_header.comments() {
        header.push_comment(co.as_bytes());
    }

    // Add Bramble PG tag
    let mut pg = HeaderRecord::new(b"PG");
    pg.push_tag(b"ID", "bramble-rs");
    pg.push_tag(b"VN", env!("CARGO_PKG_VERSION"));
    pg.push_tag(b"CL", pg_args);
    header.push_record(&pg);

    header
}
