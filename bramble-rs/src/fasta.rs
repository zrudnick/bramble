use anyhow::Result;
use needletail::parse_fastx_file;
use std::collections::HashMap;
use std::path::Path;

#[derive(Debug, Default)]
pub struct FastaDb {
    seqs: HashMap<String, Vec<u8>>,
}

impl FastaDb {
    /// Build a `FastaDb` from already-extracted reference sequences, keyed by
    /// contig name (the same first-whitespace-token convention `load` produces).
    ///
    /// Lets callers populate the rescue reference from a source other than a
    /// FASTA file on disk — e.g. the reference sequences already resident in a
    /// loaded aligner index (rammap/minimap2) — so genome-mode soft-clip rescue
    /// works without a separate FASTA. Sequence bytes should be ASCII
    /// nucleotides (`get_slice` upper-cases and N-maps as needed), exactly as if
    /// they had been read from a FASTA, so downstream exon slicing is identical.
    pub fn from_seqs(seqs: HashMap<String, Vec<u8>>) -> Self {
        Self { seqs }
    }

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

    /// `from_seqs` yields a db whose `get_slice` behaves exactly like one loaded
    /// from a FASTA (1-based half-open, upper-cased) — so an index-sourced rescue
    /// reference is interchangeable with a file-sourced one.
    #[test]
    fn from_seqs_slices_like_load() {
        let mut m = HashMap::new();
        m.insert("chr1".to_string(), b"acgtACGTNN".to_vec());
        let db = FastaDb::from_seqs(m);
        assert_eq!(db.get_slice("chr1", 1, 5).unwrap(), b"ACGT".to_vec());
        assert_eq!(db.get_slice("chr1", 5, 9).unwrap(), b"ACGT".to_vec());
        assert!(db.get_slice("missing", 1, 5).is_none());
    }

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
