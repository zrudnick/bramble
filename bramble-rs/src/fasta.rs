use anyhow::Result;
use needletail::parse_fastx_file;
use std::collections::HashMap;
use std::path::Path;

#[derive(Debug, Default)]
pub struct FastaDb {
    seqs: HashMap<String, Vec<u8>>,
}

impl FastaDb {
    pub fn load(path: &Path) -> Result<Self> {
        let mut reader = parse_fastx_file(path)
            .map_err(|e| anyhow::anyhow!("failed to open FASTA {}: {}", path.display(), e))?;
        let mut seqs: HashMap<String, Vec<u8>> = HashMap::new();

        while let Some(result) = reader.next() {
            let record = result
                .map_err(|e| anyhow::anyhow!("failed to parse FASTA record: {}", e))?;
            let name = std::str::from_utf8(record.id())
                .unwrap_or("")
                .to_string();
            let seq = record.seq().to_vec();
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
