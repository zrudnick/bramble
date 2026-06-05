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
            // needletail's `id()` returns the entire FASTA header line after `>`
            // (description included), but reference names everywhere else — GTF/GFF
            // `seqname`, BAM `@SQ`/`get_slice` lookups, samtools faidx, minimap2 —
            // use only the first whitespace-delimited token (the accession). Key by
            // that token so lookups resolve; otherwise every `get_slice` misses and
            // exon sequences come back empty, silently disabling soft-clip rescue.
            let header = std::str::from_utf8(record.id()).unwrap_or("");
            let name = header.split_whitespace().next().unwrap_or("").to_string();
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    /// Real genome FASTAs (RefSeq/Ensembl/GENCODE) carry a description after the
    /// accession in the header. The accession is the name used everywhere else
    /// (GTF `seqname`, BAM `@SQ`, `get_slice` calls), so we must key on the first
    /// whitespace-delimited token. Regression guard for the bug where keying on
    /// the full header silently emptied every exon sequence and disabled
    /// soft-clip rescue.
    #[test]
    fn keys_by_first_header_token_not_full_description() {
        let dir = std::env::temp_dir();
        let path = dir.join(format!("bramble_fasta_keytest_{}.fa", std::process::id()));
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, ">NC_000001.11 Homo sapiens chromosome 1, GRCh38.p14").unwrap();
            writeln!(f, "ACGTACGTAC").unwrap();
            writeln!(f, ">chr2\tsecond contig with tab").unwrap();
            writeln!(f, "TTTTGGGGCC").unwrap();
        }
        let db = FastaDb::load(&path).unwrap();
        // Lookup by the bare accession (as a GTF seqname would supply) must hit.
        assert_eq!(db.get_slice("NC_000001.11", 1, 5).unwrap(), b"ACGT".to_vec());
        // Tab-delimited descriptions are stripped too.
        assert_eq!(db.get_slice("chr2", 1, 5).unwrap(), b"TTTT".to_vec());
        // The full header line must NOT be a valid key.
        assert!(db.get_slice("NC_000001.11 Homo sapiens chromosome 1, GRCh38.p14", 1, 5).is_none());
        std::fs::remove_file(&path).ok();
    }
}
