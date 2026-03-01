// bam_input.rs is used only by the binary (BAM I/O path).
#![allow(dead_code)]
use crate::types::{HashMap, RefId};
use anyhow::Result;
use rust_htslib::bam;
use rust_htslib::bam::Read as HtsRead;
use std::path::Path;

pub struct BamInput {
    pub refname_to_id: HashMap<String, RefId>,
    pub reader: bam::Reader,
}

pub fn open_bam(path: &Path) -> Result<BamInput> {
    let reader = bam::Reader::from_path(path)?;
    let refname_to_id = {
        let header = reader.header();
        header
            .target_names()
            .iter()
            .enumerate()
            .map(|(i, n)| (String::from_utf8_lossy(n).to_string(), i as RefId))
            .collect()
    };
    Ok(BamInput { refname_to_id, reader })
}
