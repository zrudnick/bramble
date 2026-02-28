use anyhow::{anyhow, Result};
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InputFormat {
    Gtf,
    Gff3,
}

#[derive(Debug, Clone)]
pub struct Exon {
    pub start: u32,
    pub end: u32,
}

#[derive(Debug, Clone)]
pub struct Transcript {
    pub id: String,
    pub seqname: String,
    pub strand: char,
    pub exons: Vec<Exon>,
}

impl Transcript {
    pub fn length(&self) -> u32 {
        self.exons
            .iter()
            .map(|e| e.end.saturating_sub(e.start))
            .sum()
    }
}

pub fn detect_format(path: &Path) -> Result<InputFormat> {
    let ext = path
        .extension()
        .and_then(|s| s.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();
    match ext.as_str() {
        "gtf" => Ok(InputFormat::Gtf),
        "gff" | "gff3" => Ok(InputFormat::Gff3),
        _ => Err(anyhow!(
            "unable to detect annotation format from extension: .{}",
            ext
        )),
    }
}

/// Load transcript and exon features from GTF/GFF.
///
/// Coordinate conventions:
/// - GTF/GFF are 1-based inclusive.
/// - To match the original C++ implementation (gclib/GSamRecord),
///   we store genomic intervals as 1-based, half-open [start, end+1).
///   That means we keep `start` as-is and convert `end` -> `end + 1`.
pub fn load_transcripts(path: &Path) -> Result<Vec<Transcript>> {
    match detect_format(path)? {
        InputFormat::Gtf => load_gtf(path),
        InputFormat::Gff3 => load_gff3(path),
    }
}

fn load_gtf(path: &Path) -> Result<Vec<Transcript>> {
    // NOTE: We parse with noodles-gtf, but we only extract fields and attributes.
    // record_bufs yields gff::feature::RecordBuf, which provides a uniform API.
    let reader = File::open(path)?;
    let mut reader = noodles::gtf::io::Reader::new(BufReader::new(reader));

    let mut transcripts: HashMap<String, Transcript> = HashMap::new();

    for result in reader.record_bufs() {
        let record = result?;

        // Only use transcript + exon features
        let feature_type: &[u8] = record.ty().as_ref();
        if feature_type != b"transcript" && feature_type != b"exon" {
            continue;
        }

        let seqname = record.reference_sequence_name().to_string();
        let strand = strand_to_char(record.strand());

        // Convert 1-based inclusive -> 0-based half-open
        let start_1 = record.start().get();
        let end_1 = record.end().get();
        let start = u32::try_from(start_1).map_err(|_| anyhow!("GTF start out of range"))?;
        let end = u32::try_from(end_1.saturating_add(1))
            .map_err(|_| anyhow!("GTF end out of range"))?;

        let attrs = record.attributes();
        let transcript_id = get_record_buf_attribute(attrs, b"transcript_id")
            .ok_or_else(|| anyhow!("missing transcript_id in GTF attributes"))?;

        let entry = transcripts.entry(transcript_id.clone()).or_insert_with(|| Transcript {
            id: transcript_id.clone(),
            seqname: seqname.clone(),
            strand,
            exons: Vec::new(),
        });

        if feature_type == b"exon" {
            entry.exons.push(Exon { start, end });
        }
    }

    Ok(transcripts.into_values().collect())
}

fn load_gff3(path: &Path) -> Result<Vec<Transcript>> {
    let reader = File::open(path)?;
    let mut reader = noodles::gff::io::Reader::new(BufReader::new(reader));

    let mut transcripts: HashMap<String, Transcript> = HashMap::new();

    for result in reader.record_bufs() {
        let record = result?;

        let feature_type: &[u8] = record.ty().as_ref();
        if feature_type != b"transcript" && feature_type != b"exon" {
            continue;
        }

        let seqname = record.reference_sequence_name().to_string();
        let strand = strand_to_char(record.strand());

        let start_1 = record.start().get();
        let end_1 = record.end().get();
        let start = u32::try_from(start_1).map_err(|_| anyhow!("GFF3 start out of range"))?;
        let end = u32::try_from(end_1.saturating_add(1))
            .map_err(|_| anyhow!("GFF3 end out of range"))?;

        let attrs = record.attributes();
        let transcript_id = if feature_type == b"transcript" {
            get_record_buf_attribute(attrs, b"ID")
        } else {
            get_record_buf_attribute(attrs, b"Parent")
        }
        .ok_or_else(|| anyhow!("missing transcript id in GFF3 attributes"))?;

        let entry = transcripts.entry(transcript_id.clone()).or_insert_with(|| Transcript {
            id: transcript_id.clone(),
            seqname: seqname.clone(),
            strand,
            exons: Vec::new(),
        });

        if feature_type == b"exon" {
            entry.exons.push(Exon { start, end });
        }
    }

    Ok(transcripts.into_values().collect())
}

fn get_record_buf_attribute(
    attrs: &noodles::gff::feature::record_buf::Attributes,
    key: &[u8],
) -> Option<String> {
    let value = attrs.get(key)?;
    value.iter().next().map(|v| v.to_string())
}

fn strand_to_char(strand: noodles::gff::feature::record::Strand) -> char {
    use noodles::gff::feature::record::Strand;
    match strand {
        Strand::Forward => '+',
        Strand::Reverse => '-',
        Strand::None => '.',
        Strand::Unknown => '?',
    }
}
