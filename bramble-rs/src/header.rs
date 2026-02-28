use crate::annotation::Transcript;
use noodles::sam::{
    self,
    header::record::value::{
        map::{self, header::Version, ReferenceSequence},
        Map,
    },
};
use std::num::NonZeroUsize;

pub fn build_header(transcripts: &[Transcript]) -> sam::Header {
    let mut header = sam::Header::builder();

    let hd = Map::<map::Header>::builder()
        .set_version(Version::new(1, 0))
        .build()
        .unwrap_or_default();
    header = header.set_header(hd);

    let mut refs = sam::header::ReferenceSequences::default();
    for tx in transcripts {
        let len = tx.length() as usize;
        let Some(nz) = NonZeroUsize::new(len) else {
            continue;
        };
        let rs = Map::<ReferenceSequence>::new(nz);
        refs.insert(tx.id.clone().into(), rs);
    }

    header = header.set_reference_sequences(refs);
    header.build()
}
