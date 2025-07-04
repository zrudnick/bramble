// reads.cpp

#include "bramble.h"
#include "htslib/sam.h"

extern bool LONG_READS;
static GStr read_id("",
                    256); // Prevent repeated reallocation for each parsed read

void process_exons(GSamRecord *brec, CReadAln *readaln) {
  auto exons = brec->exons;
  for (int i = 0; i < exons.Count(); i++) {
    readaln->len += exons[i].len();
    readaln->segs.Add(exons[i]);
  }
}

int add_new_read(BundleData &bundle, CReadAln *readaln, GSamRecord *brec) {
  int n = bundle.reads.Add(readaln);
  bundle.reads[n]->brec = brec;
  return n;
}

void update_bundle_end(uint bundle_end, BundleData &bundle, int read_end) {
  if (read_end > bundle_end) {
    bundle_end = read_end;
    bundle.end = bundle_end;
  }
}

float calculate_read_count(GSamRecord *brec) {
  float unitig_cov = brec->tag_float("YK");
  float read_count = static_cast<float>(brec->tag_int("YC"));
  if (read_count == 0)
    read_count = 1;
  if (unitig_cov)
    read_count = unitig_cov;
  return read_count;
}

std::string create_read_id(const char *read_name, int position, int hi) {
  std::string id(read_name);
  id += '-';
  id += position;
  id += ".=";
  id += hi;
  return id;
}

void add_pair_if_new(CReadAln *read_aln, int pair_index, float read_count) {
  // Check if pairing already exists
  for (int i = 0; i < read_aln->pair_idx.Count(); i++) {
    if (read_aln->pair_idx[i] == pair_index) {
      return; // Pairing already exists
    }
  }

  // Add new pairing
  read_aln->pair_idx.Add(pair_index);
  read_aln->pair_count.Add(read_count);
}

void establish_pairing(BundleData &bundle, int read1_index, int read2_index,
                       float read_count) {

  // Add pairing information to both reads if not already present
  add_pair_if_new(bundle.reads[read1_index], read2_index, read_count);
  add_pair_if_new(bundle.reads[read2_index], read1_index, read_count);
}

void process_paired_reads(BundleData &bundle, int bundle_start, int read_start,
                          int read_index, GSamRecord *brec, float read_count,
                          int hi, GHash<int> &hashread) {

  // Only process pairs on same chromosome/contig
  if (brec->refId() != brec->mate_refId())
    return;

  int mate_start = brec->mate_start();

  // Ignore pairs in previous bundles
  if (mate_start < bundle_start)
    return;

  // Handle pair processing based on mate position
  if (mate_start <= read_start) {
    std::string read_id =
        create_read_id(bundle.reads[read_index]->brec->name(), mate_start, hi);

    const int *mate_index = hashread[read_id.c_str()];
    if (mate_index) {
      establish_pairing(bundle, read_index, *mate_index, read_count);
      hashread.Remove(read_id.c_str());
    }

  } else {
    std::string read_id = create_read_id(brec->name(), read_start, hi);
    hashread.Add(read_id.c_str(), read_index);
  }
}

void process_read_in(uint bundle_start, uint bundle_end, BundleData &bundle,
                     GHash<int> &hashread, GSamRecord *brec, char strand, int nh, int hi) {

  // Skip secondary alignments
  // if (brec->flags() & BAM_FSECONDARY) return;

  // Extract alignment information
  int read_start = brec->start;

  // Create new read alignment
  CReadAln *readaln = new CReadAln(strand, nh, brec->start, brec->end);
  readaln->longread = (LONG_READS || brec->uval);

  // Process exons
  process_exons(brec, readaln);

  // Add read to bundle
  int n = add_new_read(bundle, readaln, brec);

  // Update bundle end if necessary
  update_bundle_end(bundle_end, bundle, brec->end);

  // Calculate read count with multi-mapping correction
  float read_count = calculate_read_count(brec);

  // Handle paired-end reads
  process_paired_reads(bundle, bundle_start, read_start, n, brec, read_count,
                       hi, hashread);
}
