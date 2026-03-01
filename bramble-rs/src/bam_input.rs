// bam_input.rs is used only by the binary (BAM I/O path).
#![allow(dead_code)]
use crate::types::{HashMap, HashMapExt, RefId};
use anyhow::Result;
use noodles::{bam, bgzf, sam};
use std::fs::File;
use std::path::Path;

pub struct BamInput {
    pub refname_to_id: HashMap<String, RefId>,
    pub reader: bam::io::Reader<bgzf::io::Reader<File>>,
}

pub fn open_bam(path: &Path) -> Result<BamInput> {
    let file = File::open(path)?;
    let mut reader = bam::io::Reader::new(file);

    let header = reader.read_header()?;
    let refname_to_id = build_refid_map(header.reference_sequences());

    Ok(BamInput { refname_to_id, reader })
}

fn build_refid_map(
    reference_sequences: &sam::header::ReferenceSequences,
) -> HashMap<String, RefId> {
    let mut map = HashMap::new();
    for (i, (name, _rs)) in reference_sequences.iter().enumerate() {
        map.insert(name.to_string(), i as RefId);
    }
    map
}
