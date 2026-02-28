use anyhow::Result;
use noodles::fasta;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

#[derive(Debug, Default)]
pub struct FastaDb {
    seqs: HashMap<String, Vec<u8>>,
}

impl FastaDb {
    pub fn load(path: &Path) -> Result<Self> {
        let file = File::open(path)?;
        let mut reader = fasta::io::Reader::new(BufReader::new(file));
        let mut seqs: HashMap<String, Vec<u8>> = HashMap::new();

        for result in reader.records() {
            let record = result?;
            let name = String::from_utf8_lossy(record.name()).to_string();
            let seq = record.sequence().as_ref().to_vec();
            seqs.insert(name, seq);
        }

        Ok(Self { seqs })
    }

    pub fn get_slice(&self, seqname: &str, start: u32, end: u32) -> Option<Vec<u8>> {
        let seq = self.seqs.get(seqname)?;
        if start == 0 || end == 0 {
            return None;
        }
        // Exon coordinates are stored as 1-based, half-open [start, end),
        // but Rust slices are 0-based, half-open. Match the C++ copyRange(start, end-1).
        let s = start.saturating_sub(1) as usize;
        let e = end.saturating_sub(1) as usize;
        if s <= e && e <= seq.len() {
            let mut out = seq[s..e].to_vec();
            for b in &mut out {
                *b = b.to_ascii_uppercase();
            }
            Some(out)
        } else {
            None
        }
    }
}
