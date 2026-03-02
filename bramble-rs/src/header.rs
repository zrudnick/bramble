// header.rs is used only by the binary.
#![allow(dead_code)]
use crate::annotation::Transcript;

pub fn build_hts_header(transcripts: &[Transcript]) -> rust_htslib::bam::Header {
    use rust_htslib::bam::header::HeaderRecord;
    let mut header = rust_htslib::bam::Header::new();
    let mut hd = HeaderRecord::new(b"HD");
    hd.push_tag(b"VN", "1.0");
    header.push_record(&hd);
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
    header
}
