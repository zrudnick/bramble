// header.rs is used only by the binary.
#![allow(dead_code)]
use crate::annotation::Transcript;
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

pub fn build_hts_header(transcripts: &[Transcript], input_header: &Header, pg_args: &str) -> Header {
    let mut header = Header::new();

    // HD record
    let mut hd = HeaderRecord::new(b"HD");
    hd.push_tag(b"VN", "1.0");
    hd.push_tag(b"SO", "coordinate");
    header.push_record(&hd);

    // SQ records from transcripts
    for tx in transcripts {
        let len = tx.length();
        if len == 0 { continue; }
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", &tx.id);
        sq.push_tag(b"LN", len);
        header.push_record(&sq);
    }

    // Retain all PG records from input BAM header
    let header_view = input_header.to_hashmap();
    if let Some(pg_records) = header_view.get("PG") {
        for pg in pg_records {
            let mut rec = HeaderRecord::new(b"PG");
            for (key, val) in pg {
                rec.push_tag(key.as_bytes(), val);
            }
            header.push_record(&rec);
        }
    }

    // Add Bramble PG tag
    let mut pg = HeaderRecord::new(b"PG");
    pg.push_tag(b"ID", "bramble-rs");
    pg.push_tag(b"VN", env!("CARGO_PKG_VERSION"));
    pg.push_tag(b"CL", pg_args);
    header.push_record(&pg);

    header
}
