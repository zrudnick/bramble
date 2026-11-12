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

/// Program name in our @PG lines (ID and PN).
const PG_NAME: &str = "bramble-rs";

/// The @PG lines to add, as htslib's `sam_hdr_add_pg` makes them (the C++
/// version calls it): one line per chain end of the input's @PG lines, with PP
/// set to that end, or a single line without PP when there are none. ID is the
/// program name, or `name.N` (N counting from 1) when already in use. Tag order
/// as htslib: ID, PN, PP, then VN and CL.
fn new_pg_records(
    input_pg_ids: &[(String, Option<String>)],
    pg_args: &str,
) -> Vec<HeaderRecord<'static>> {
    let mut used: Vec<String> = input_pg_ids.iter().map(|(id, _)| id.clone()).collect();
    // chain ends: @PG lines no other @PG line names as its PP, in header order
    let ends: Vec<String> = input_pg_ids
        .iter()
        .filter(|(id, _)| {
            !input_pg_ids
                .iter()
                .any(|(_, pp)| pp.as_deref() == Some(id.as_str()))
        })
        .map(|(id, _)| id.clone())
        .collect();
    let mut counter = 1;
    let mut next_id = |used: &mut Vec<String>| {
        let mut id = PG_NAME.to_string();
        while used.contains(&id) {
            id = format!("{PG_NAME}.{counter}");
            counter += 1;
        }
        used.push(id.clone());
        id
    };
    let make = |id: String, pp: Option<&str>| {
        let mut pg = HeaderRecord::new(b"PG");
        pg.push_tag(b"ID", id);
        pg.push_tag(b"PN", PG_NAME);
        if let Some(pp) = pp {
            pg.push_tag(b"PP", pp.to_string());
        }
        pg.push_tag(b"VN", env!("CARGO_PKG_VERSION"));
        pg.push_tag(b"CL", pg_args.to_string());
        pg
    };
    if ends.is_empty() {
        return vec![make(next_id(&mut used), None)];
    }
    ends.iter()
        .map(|end| {
            let id = next_id(&mut used);
            make(id, Some(end))
        })
        .collect()
}

pub fn build_hts_header(
    transcripts: &[Transcript],
    input_header: &Header,
    pg_args: &str,
    guide_path: &str,
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
    let input_pg_ids: Vec<(String, Option<String>)> = header_view
        .get("PG")
        .map(|records| {
            records
                .iter()
                .filter_map(|r| r.get("ID").map(|id| (id.clone(), r.get("PP").cloned())))
                .collect()
        })
        .unwrap_or_default();
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

    // Our @PG line(s), chained to the input's
    for pg in new_pg_records(&input_pg_ids, pg_args) {
        header.push_record(&pg);
    }
    // the annotation the transcripts came from (C++ writes the same line)
    header.push_comment(format!("Generated from GTF: {guide_path}").as_bytes());

    header
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pg_lines(ids: &[(&str, Option<&str>)]) -> Vec<String> {
        let ids: Vec<(String, Option<String>)> = ids
            .iter()
            .map(|(id, pp)| (id.to_string(), pp.map(|p| p.to_string())))
            .collect();
        new_pg_records(&ids, "bramble-rs in.bam")
            .iter()
            .map(|r| {
                let mut h = Header::new();
                h.push_record(r);
                String::from_utf8(h.to_bytes())
                    .unwrap()
                    .trim_end()
                    .to_string()
            })
            .collect()
    }

    const VN: &str = env!("CARGO_PKG_VERSION");

    #[test]
    fn pg_without_input_programs_has_no_pp() {
        assert_eq!(
            pg_lines(&[]),
            vec![format!(
                "@PG\tID:bramble-rs\tPN:bramble-rs\tVN:{VN}\tCL:bramble-rs in.bam"
            )]
        );
    }

    #[test]
    fn pg_chains_to_the_last_program() {
        let lines = pg_lines(&[("hisat2", None), ("samtools", Some("hisat2"))]);
        assert_eq!(lines.len(), 1);
        assert!(lines[0].starts_with("@PG\tID:bramble-rs\tPN:bramble-rs\tPP:samtools\t"));
    }

    #[test]
    fn pg_id_made_unique_and_one_line_per_chain() {
        // two chains: hisat2 -> bramble-rs, and minimap2
        let lines = pg_lines(&[
            ("hisat2", None),
            ("bramble-rs", Some("hisat2")),
            ("minimap2", None),
        ]);
        assert_eq!(lines.len(), 2);
        assert!(lines[0].starts_with("@PG\tID:bramble-rs.1\tPN:bramble-rs\tPP:bramble-rs\t"));
        assert!(lines[1].starts_with("@PG\tID:bramble-rs.2\tPN:bramble-rs\tPP:minimap2\t"));
    }
}
