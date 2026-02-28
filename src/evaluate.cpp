
#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include "types.h"
#include "bramble.h"
#include "g2t.h"
#include "evaluate.h"
#include "bam.h"
#include "ksw2.h"

#ifndef NOTHREADS
#include "GThreads.h"
#endif

extern bool BRAMBLE_DEBUG;
extern bool LONG_READS;
extern bool USE_FASTA;
extern bool SOFT_CLIPS;
extern bool STRICT;

namespace bramble {

  ReadEvaluator::~ReadEvaluator() {}

  // Evaluate read group and return matches
  std::vector<ExonChainMatch> 
  ReadEvaluator::evaluate(CReadAln &read,
                          std::shared_ptr<g2tTree> g2t,
                          uint8_t *seq, int seq_len) {
    return std::vector<ExonChainMatch>{};
  }

  ExonStatus get_exon_status(uint32_t exon_count, uint32_t j) {
    ExonStatus status;
    if (exon_count == 1) {
      status = ONLY_EXON;
    } else if (j == 0) {
      status = FIRST_EXON;
    } else if (j < (exon_count - 1)) {
      status = MIDDLE_EXON;
    } else {
      status = LAST_EXON;  
    }
    return status;
  }

  std::vector<char> 
  ReadEvaluator::get_strands_to_check(CReadAln &read) {
    // read strand is an informed guess as to which
    // strand the transcript could have come from
    if (LONG_READS) return {'+', '-'};

    if (read.strand == '+') return {'+'}; 
    if (read.strand == '-') return {'-'}; 
    return {'+', '-'};  // unstranded: try forward first, then reverse
  }

  void ReadEvaluator::get_clips(CReadAln &read, ReadEvaluationConfig &config,
                                bool &failure, bool &has_left_clip, bool &has_right_clip,
                                uint32_t &n_left_clip, uint32_t &n_right_clip) {
    bam1_t* b = read.brec->get_b();
    uint32_t* cigar = bam_get_cigar(b);
    uint32_t n_cigar = b->core.n_cigar;

    if (!cigar) {
      printf("Error: alignment is missing CIGAR string\n");
      printf("Skipping this alignment\n");
      failure = true;
      return;
    }
    
    if ((cigar[0] & BAM_CIGAR_MASK) == BAM_CHARD_CLIP) {
      if (1 < n_cigar && ((cigar[1] & BAM_CIGAR_MASK) == BAM_CSOFT_CLIP)) {
        has_left_clip = USE_FASTA;
        n_left_clip = cigar[1] >> BAM_CIGAR_SHIFT;
      }
       
    } else if ((cigar[0] & BAM_CIGAR_MASK) == BAM_CSOFT_CLIP) {
      has_left_clip = USE_FASTA;
      n_left_clip = cigar[0] >> BAM_CIGAR_SHIFT;
    }
      
    if ((cigar[n_cigar - 1] & BAM_CIGAR_MASK) == BAM_CHARD_CLIP) {
      if (n_cigar - 2 >= 0 && ((cigar[n_cigar - 2] & BAM_CIGAR_MASK) == BAM_CSOFT_CLIP)) {
        has_right_clip = USE_FASTA;
        n_right_clip = cigar[n_cigar - 2] >> BAM_CIGAR_SHIFT;
      }
        
    } else if ((cigar[n_cigar - 1] & BAM_CIGAR_MASK) == BAM_CSOFT_CLIP) {
      has_right_clip = USE_FASTA;
      n_right_clip = (cigar[n_cigar - 1] >> BAM_CIGAR_SHIFT);
    }

    if (config.print) {
      printf("Has left clip? %d # left clip = %d\n", has_left_clip, n_left_clip);
      printf("Has right clip? %d # right clip = %d\n", has_right_clip, n_right_clip);
    }
  }

  bool ReadEvaluator::correct_for_gaps(TidData &td, tid_t tid,
                                      Segment &seg2,
                                      ReadEvaluationConfig &config, 
                                      std::shared_ptr<g2tTree> g2t,
                                      char strand, refid_t refid) {
    size_t n_segs = td.segments.size();
    const auto &seg1 = td.segments[n_segs - 1];

    uint8_t gap = seg2.gexon.exon_id - seg1.gexon.exon_id;

    // short reads: require strict continuity
    if (!LONG_READS) {
      if (gap != 1) {
        td.elim = true;
        return false;
      } else {
        return true;
      }
    }

    else { // long reads: handle gaps
      if (gap > 2) {
        td.elim = true;
        return false;
      }

      else if (gap == 2) {
        uint32_t gap_start = (strand == '+') ? 
          seg2.gexon.prev_start : seg2.gexon.next_start;
        uint32_t gap_end = (strand == '+') ? 
          seg2.gexon.prev_end : seg2.gexon.next_end;

        if ((gap_start == 0 && gap_end == 0) || 
          (gap_end - gap_start > config.max_gap)) {
          td.elim = true;
          return false;
        }

        GuideExon gap_exon;
        if (!g2t->getGuideExonForTid(refid, 
          strand, tid, gap_start, gap_end, gap_exon)) {
          td.elim = true;
          return false;
        }

        Segment gap_seg;
        gap_seg.gexon = gap_exon;
        gap_seg.has_gexon = true;
        gap_seg.status = GAP_EXON;
        gap_seg.is_small_exon = (gap_exon.end - gap_exon.start 
          <= config.small_exon_size);
        
        td.segments.emplace_back(gap_seg);
      }

      return true;
    }
  }

  void ReadEvaluator::get_intervals(unordered_map<tid_t, TidData> &data,
                                    CReadAln &read, uint32_t j, uint32_t exon_count,
                                    ReadEvaluationConfig &config, std::shared_ptr<g2tTree> g2t,
                                    refid_t refid, char strand, bool has_left_clip,
                                    bool has_right_clip, bool &failure) {
    const auto &qexon = read.segs[j];

    auto status = get_exon_status(exon_count, j);

    bool is_small_exon = (qexon.end - qexon.start <= config.small_exon_size);
    bool data_empty = data.empty();

    thread_local std::vector<GuideExon> guide_exons;
    guide_exons.clear();

    std::vector<tid_t> candidate_tids;
    candidate_tids.reserve(100);

    if (g2t->getGuideExons(refid, strand, 
      qexon, config, status, guide_exons)) {

      for (auto &gexon : guide_exons) {
        tid_t tid = gexon.tid;
        candidate_tids.emplace_back(tid);

        Segment segment{ 
          .has_gexon = true,
          .has_qexon = true,
          .gexon = gexon,
          .qexon = qexon, 
          .status = status, 
          .is_small_exon = is_small_exon 
        };

        if (data_empty) {
          TidData td;
          td.has_left_clip = has_left_clip;
          td.has_right_clip = has_right_clip;
          td.segments.emplace_back(segment);

          data[tid] = td;
        } else {

          auto it = data.find(tid); 
          if (it == data.end() || it->second.elim) {
            continue;
          }
          TidData& td = it->second;
          
          correct_for_gaps(td, tid, segment, config, g2t, strand, refid);

          td.segments.emplace_back(segment);
        }
      }

      for (auto &pair : data) {
        tid_t tid = pair.first;
        if (std::find(candidate_tids.begin(), 
          candidate_tids.end(), tid) == candidate_tids.end()) {
          pair.second.elim = true;
        }
      }

      return;
    }

    if (status != ONLY_EXON && config.ignore_small_exons 
      && is_small_exon) {
    
      if (status == MIDDLE_EXON) {
        if (data.empty()) {
          failure = true;
          return;
        }

        for (auto &pair: data) {
          tid_t tid = pair.first;
          auto &td = pair.second;
          
          Segment ignore;
          ignore.qexon = qexon;
          ignore.has_qexon = true;
          ignore.has_gexon = false;
          ignore.status = INS_EXON;
          ignore.is_small_exon = true;
          td.segments.emplace_back(ignore);
        }

        // success
        return;
      } else {
        failure = true;
        return;
      }
    }

    failure = true;
    return;
  }

  // from https://github.com/lh3/ksw2/tree/master
  static kswResult align(const char *tseq, const char *qseq, int sc_mch, 
                        int sc_mis, int gapo, int gape, int zdrop) {
    kswResult result;
    int i, a = sc_mch, b = sc_mis < 0? sc_mis : -sc_mis;
    int8_t mat[25] = { (int8_t)a,(int8_t)b,(int8_t)b,(int8_t)b,0,
                       (int8_t)b,(int8_t)a,(int8_t)b,(int8_t)b,0,
                       (int8_t)b,(int8_t)b,(int8_t)a,(int8_t)b,0,
                       (int8_t)b,(int8_t)b,(int8_t)b,(int8_t)a,0,
                       0,0,0,0,0 };
    int tl = strlen(tseq), ql = strlen(qseq);
    uint8_t *ts, *qs, c[256];
    ksw_extz_t ez;

    memset(&ez, 0, sizeof(ksw_extz_t));
    memset(c, 4, 256);
    c['A'] = c['a'] = 0; c['C'] = c['c'] = 1;
    c['G'] = c['g'] = 2; c['T'] = c['t'] = 3;
    ts = (uint8_t*)malloc(tl);
    qs = (uint8_t*)malloc(ql);
    for (i = 0; i < tl; ++i) ts[i] = c[(uint8_t)tseq[i]];
    for (i = 0; i < ql; ++i) qs[i] = c[(uint8_t)qseq[i]];
    
    int flag = KSW_EZ_EXTZ_ONLY & KSW_EZ_APPROX_MAX & KSW_EZ_APPROX_DROP;
    ksw_extz2_sse(0, ql, qs, tl, ts, 5, mat, gapo, gape, -1, zdrop, 0, flag, &ez);
    
    result.score = ez.score;
    result.max = ez.max;  // use max due to EXTZ_ONLY
    result.cigar_array = ez.cigar;
    result.n_cigar = ez.n_cigar;
    
    free(ts); free(qs);
    return result;
  }

  // Extracts the query sequence for left-clip rescue
  static std::string 
  get_left_rescue_query(uint8_t *seq, uint32_t n_left_clip, 
                        uint32_t left_ins) {
    uint32_t total = n_left_clip + left_ins;
    std::string q;
    q.reserve(total);
    for (uint32_t i = 0; i < total; i++)
      q += seq_nt16_str[bam_seqi(seq, i)];
    return q;
  }

  // Concatenate guide exon sequences
  static bool 
  collect_left_exons(const std::string &qseq, std::shared_ptr<g2tTree> g2t,
                    refid_t refid, char strand, tid_t tid, const GuideExon &start,
                    std::string &out_gseq, std::vector<GuideExon> &out_exons, bool print) {
    GuideExon curr = start;
    int i = 0;

    while (qseq.length() > out_gseq.length()) {
      i++;

      bool has_neighbor = (strand == '+') ? curr.has_prev : curr.has_next;
      if (!has_neighbor) {
        if (i == 1) {
          if (print) printf("no neighbor on left side :(\n");
          return false;
        }
        break;
      }

      GuideExon next;
      if (strand == '+')
        g2t->getGuideExonForTid(refid, strand, tid, curr.prev_start, curr.prev_end, next);
      else
        g2t->getGuideExonForTid(refid, strand, tid, curr.next_start, curr.next_end, next);

      std::string next_seq = next.seq;   
      out_gseq = next_seq + out_gseq;
      out_exons.insert(out_exons.begin(), next);
      curr = next;
    }

    return !out_exons.empty();
  }

  // from https://github.com/lh3/ksw2/tree/master
  static kswResult 
  align_reversed(const std::string &qseq, const std::string &gseq, bool print) {
    int len1 = qseq.length();
    int len2 = gseq.length();

    // get a window 40bp larger than qseq
    int start_pos = std::max(0, len2 - (len1 + 40));
    std::string gseq_short = gseq.substr(start_pos);

    // reverse both so ksw2 extension runs right-to-left
    std::string qrev = qseq, grev = gseq_short;
    std::reverse(qrev.begin(), qrev.end());
    std::reverse(grev.begin(), grev.end());

    if (print) {
      printf("\nRunning alignment:\n");
      printf("qseq = %s (len=%zu)\n", qrev.c_str(), qrev.length());
      printf("gseq = %s (len=%zu)\n", grev.c_str(), grev.length());
    }

    return align(grev.c_str(), qrev.c_str(),
      1,    // match score
      -4,   // mismatch penalty
      4,    // gap open
      1,    // gap extend
      40   // z-drop
    );
  }

  static Segment build_left_clip_segment(const kswResult &result, int q_len, int small_exon) {
    Segment seg = {};
    seg.has_qexon     = false;
    seg.has_gexon     = false;
    seg.status        = LEFTC_EXON;
    seg.is_small_exon = (q_len <= small_exon);
    seg.cigar         = Cigar();
    seg.score         = result.max;  // use max, not score, for extension alignment

    // CIGAR was built for the reversed sequences, so iterate in reverse
    for (int i = result.n_cigar - 1; i >= 0; --i) {
      int  len     = bam_cigar_oplen(result.cigar_array[i]);
      int  op      = bam_cigar_op(result.cigar_array[i]);
      char op_char = "MID"[op];

      if      (i == result.n_cigar - 1 && op_char == 'D') { /* leading deletion — discard */ }
      else if (i == result.n_cigar - 1 && op_char == 'I') seg.cigar.add_operation(len, BAM_CLIP_OVERRIDE);
      else if (op_char == 'D')                            seg.cigar.add_operation(len, BAM_CDEL_OVERRIDE);
      else if (op_char == 'I')                            seg.cigar.add_operation(len, BAM_CINS_OVERRIDE);
      else                                                seg.cigar.add_operation(len, BAM_CMATCH_OVERRIDE);
    }

    return seg;
  }


  void ReadEvaluator::left_clip_rescue(TidData &td, char strand, 
                                      std::shared_ptr<g2tTree> g2t,
                                      refid_t refid, tid_t tid, uint32_t n_left_clip,
                                      ReadEvaluationConfig &config, CReadAln &read, 
                                      uint8_t *seq, int seq_len) {
    td.has_left_clip = false;
    Segment &seg = td.segments[0];

    if (!seg.has_gexon || seg.gexon.left_gap > 0) {
      if (seg.has_gexon && config.print)
        printf("left gap = %d\n", seg.gexon.left_gap);
      return;
    }
    auto &gexon = seg.gexon;

    std::string qseq = get_left_rescue_query(seq, n_left_clip, gexon.left_ins);

    std::string gseq;
    std::vector<GuideExon> collected_exons;
    if (!collect_left_exons(qseq, g2t, refid, strand, tid, gexon,
                            gseq, collected_exons, config.print))
      return;

    if (config.print) {
      printf("\nqseq = %s (len=%zu)\n", qseq.c_str(), qseq.length());
      printf("gseq = %s (len=%zu)\n", gseq.c_str(), gseq.length());
    }

    kswResult result = align_reversed(qseq, gseq, config.print);

    if (config.print)
      printf("\nLEFT CLIP ALIGNMENT RESULT:\nMax score = %d\n", result.max);

    if (result.max < 10 || result.score == KSW_NEG_INF) {
      free(result.cigar_array);
      return;
    }

    if (gexon.left_ins > 0) 
      gexon.left_ins = 0;

    Segment left_clip = build_left_clip_segment(result, qseq.length(), 
      config.small_exon_size);
    free(result.cigar_array);

    td.segments.emplace(td.segments.begin(), left_clip);
    td.has_left_clip = true;
  }

  static std::string 
  get_right_rescue_query(uint8_t *seq, int seq_len,
                        uint32_t n_right_clip, uint32_t right_ins) {
    uint32_t total = right_ins + n_right_clip;
    int start = seq_len - (int)total;
    std::string q;
    q.reserve(total);
    for (uint32_t i = 0; i < total; i++)
      q += seq_nt16_str[bam_seqi(seq, start + i)];
    return q;
  }

  static bool 
  collect_right_exons(const std::string &qseq, std::shared_ptr<g2tTree> g2t,
                      refid_t refid, char strand, tid_t tid,
                      const GuideExon &start, std::string &out_gseq,
                      std::vector<GuideExon> &out_exons, bool print) {
    GuideExon curr = start;
    int i = 0;

    while (qseq.length() > out_gseq.length()) {
      i++;

      bool has_neighbor = (strand == '+') ? curr.has_next : curr.has_prev;
      if (!has_neighbor) {
        if (print) printf("no neighbor on right side :(\n");
        if (i == 1) return false;
        break;
      }

      GuideExon next;
      if (strand == '+')
        g2t->getGuideExonForTid(refid, strand, tid, curr.next_start, curr.next_end, next);
      else
        g2t->getGuideExonForTid(refid, strand, tid, curr.prev_start, curr.prev_end, next);

      std::string next_seq = next.seq;
      out_gseq += next.seq;  // append (moving rightward)
      out_exons.emplace_back(next);
      curr = next;
    }

    return !out_exons.empty();
  }

  static Segment build_right_clip_segment(const kswResult &result, 
                                          int query_len, int small_exon) {
    Segment seg      = {};
    seg.has_qexon    = false;
    seg.has_gexon    = false;
    seg.status       = RIGHTC_EXON;
    seg.is_small_exon = (query_len <= small_exon);
    seg.cigar        = Cigar();
    seg.score        = result.max;  // use max, not score, for extension alignment

    for (int i = 0; i < result.n_cigar; ++i) {
      int  len     = bam_cigar_oplen(result.cigar_array[i]);
      int  op      = bam_cigar_op(result.cigar_array[i]);
      char op_char = "MID"[op];

      if      (i == (result.n_cigar - 1) && op_char == 'D') { /* trailing deletion — discard */ }
      else if (i == (result.n_cigar - 1) && op_char == 'I') seg.cigar.add_operation(len, BAM_CLIP_OVERRIDE);
      else if (op_char == 'D')                              seg.cigar.add_operation(len, BAM_CDEL_OVERRIDE);
      else if (op_char == 'I')                              seg.cigar.add_operation(len, BAM_CINS_OVERRIDE);
      else                                                  seg.cigar.add_operation(len, BAM_CMATCH_OVERRIDE);
    }

    return seg;
  }

  void ReadEvaluator::right_clip_rescue(TidData &td, char strand, 
                                        std::shared_ptr<g2tTree> g2t,
                                        refid_t refid, tid_t tid, uint32_t n_right_clip,
                                        ReadEvaluationConfig &config,
                                        CReadAln &read, uint8_t *seq, int seq_len) {
    td.has_right_clip = false;
    Segment &seg = *(td.segments.rbegin());

    if (!seg.has_gexon || seg.gexon.right_gap > 0) {
      if (seg.has_gexon && config.print)
        printf("right gap = %d\n", seg.gexon.right_gap);
      return;
    }
    auto &gexon = seg.gexon;

    std::string qseq = get_right_rescue_query(seq, seq_len, n_right_clip, gexon.right_ins);

    std::string gseq;
    std::vector<GuideExon> collected_exons;
    if (!collect_right_exons(qseq, g2t, refid, strand, tid, gexon,
                             gseq, collected_exons, config.print))
      return;

    if (config.print) {
      printf("\nRunning alignment:\n");
      printf("qseq = %s (len=%zu)\n", qseq.c_str(), qseq.length());
      printf("gseq = %s (len=%zu)\n", gseq.c_str(), gseq.length());
    }

    // Trim reference window to just past query length
    std::string gseq_short = gseq.substr(0, qseq.length() + 40);

    kswResult result = align(gseq_short.c_str(), qseq.c_str(),
      1,   // match score
      -4,  // mismatch penalty
      4,   // gap open
      1,   // gap extend
      40   // z-drop
    );

    if (config.print)
        printf("\nRIGHT CLIP ALIGNMENT RESULT:\nMax score = %d\n", result.max);

    if (result.max < 10 || result.score == KSW_NEG_INF) {
      free(result.cigar_array);
      return;
    }

    if (gexon.right_ins > 0) gexon.right_ins = 0;

    Segment right_clip = build_right_clip_segment(result, qseq.length(), config.small_exon_size);
    free(result.cigar_array);

    td.segments.emplace_back(right_clip);
    td.has_right_clip = true;
  }

  void ReadEvaluator::create_match(TidData &td, GuideExon &gexon,
                                  tid_t tid, char strand) {
    ExonChainMatch &match  = td.match;
    match.tid              = tid;
    match.align.fwpos      = gexon.pos;
    match.align.rcpos      = gexon.pos;
    match.transcript_len   = gexon.transcript_len;
    match.align.strand     = strand;
    match.align.cigar      = Cigar();
    match.align.similarity_score = 0;
    match.total_coverage   = 0;
    match.total_operations = 0;
    match.ref_consumed     = 0;
    match.prev_op          = BAM_CMATCH;
    match.junc_hits        = 0;
  }

  void ReadEvaluator::build_cigar_match(Segment &seg, TidData &td,
                                        ExonChainMatch &match, bool first_match,
                                        bool last_match, ReadEvaluationConfig &config) {
    uint32_t qstart = seg.qexon.start;
    uint32_t qend = seg.qexon.end;
    uint32_t gstart = seg.gexon.start;
    uint32_t gend = seg.gexon.end;

    uint32_t left_ins = seg.gexon.left_ins;
    uint32_t left_gap = seg.gexon.left_gap;
    uint32_t right_ins = seg.gexon.right_ins;
    uint32_t right_gap = seg.gexon.right_gap;

    if (config.print) {
      printf("qstart = %d, qend = %d\n", qstart, qend);
      printf("gstart = %d, gend = %d\n", gstart, gend);
      printf("left_ins = %d, left_gap = %d\n", left_ins, left_gap);
      printf("right_ins = %d, right_gap = %d\n", right_ins, right_gap);
    }

    Cigar &cigar = match.align.cigar;

    // Handle start boundary
    if (left_ins > 0) {
      if (seg.status == FIRST_EXON || seg.status == ONLY_EXON) {
        if (!td.has_left_clip) {
          cigar.add_operation(left_ins, BAM_CSOFT_CLIP);
          match.total_operations += left_ins;
          match.prev_op = BAM_CSOFT_CLIP;
          if (config.print) printf("Start of exon: adding softclip len %d\n", left_ins);
        }
        
      } else if (seg.status == MIDDLE_EXON || seg.status == LAST_EXON
        || td.has_left_clip) {
        cigar.add_operation(left_ins, BAM_CINS);
        match.total_operations += left_ins;
        if (match.prev_op == BAM_CDEL) {
          match.total_coverage += left_ins;
        } else if (match.prev_op == BAM_CINS) {
          match.total_operations += (match.total_operations * 0.2); // penalty for two I in a row
        }
        match.prev_op = BAM_CINS;
        if (config.print) printf("Start of exon: adding insertion len %d\n", left_ins);
      }
    } else if (left_gap > 0) {
      if (!first_match && (seg.status == MIDDLE_EXON || seg.status == LAST_EXON
        || td.has_left_clip)) {
        cigar.add_operation(left_gap, BAM_CDEL);
        match.total_operations += left_gap;
        match.ref_consumed += left_gap;
        if (match.prev_op == BAM_CINS) {
          match.total_coverage += left_gap;
        } else if (match.prev_op == BAM_CDEL) {
          match.total_operations += (match.total_operations * 0.2); // penalty for two D in a row
        }
        match.prev_op = BAM_CDEL;
        if (config.print) printf("Start of exon: adding deletion len %d\n", left_gap);
      }
    } else {
      match.junc_hits++;
    }

    // Add match for overlapping region
    uint32_t overlap_start = std::max(qstart, gstart);
    uint32_t overlap_end = std::min(qend, gend);  
    if (overlap_end >= overlap_start) {
      uint32_t match_length = overlap_end - overlap_start;
      cigar.add_operation(match_length, BAM_CMATCH);
      match.total_operations += match_length;
      match.total_coverage += match_length;
      match.ref_consumed += match_length;
      match.prev_op = BAM_CMATCH;
      if (config.print) printf("Adding match len %d\n", match_length);
    }

    // Handle end boundary
    if (right_ins > 0) {
      if (seg.status == LAST_EXON || seg.status == ONLY_EXON) {
        if (!td.has_right_clip) {
          cigar.add_operation(right_ins, BAM_CSOFT_CLIP);
          match.total_operations += right_ins;
          match.prev_op = BAM_CSOFT_CLIP;
          if (config.print) printf("End of exon: adding softclip len %d\n", right_ins);
        }
        
      } else if (seg.status == FIRST_EXON || seg.status == MIDDLE_EXON
        || td.has_right_clip) {
        cigar.add_operation(right_ins, BAM_CINS); 
        match.total_operations += right_ins;
        if (match.prev_op == BAM_CDEL) {
          match.total_coverage += right_ins;
        }
        match.prev_op = BAM_CINS;
        if (config.print) printf("End of exon: adding insertion len %d\n", right_ins);
      }
    } else if (right_gap > 0) {
      if (!last_match && (seg.status == FIRST_EXON || seg.status == MIDDLE_EXON
        || td.has_right_clip)) {
        cigar.add_operation(right_gap, BAM_CDEL); 
        match.total_operations += right_gap;
        match.ref_consumed += right_gap;
        if (match.prev_op == BAM_CINS) {
          match.total_coverage += right_gap;
        }
        match.prev_op = BAM_CDEL;
        if (config.print) printf("End of exon: adding deletion len %d\n", right_gap);
      }
    } else {
      match.junc_hits++;
    }
  }

  void ReadEvaluator::build_cigar_ins(Segment &seg, uint32_t k, uint32_t n,
                                      ExonChainMatch &match, ReadEvaluationConfig &config) {
    uint32_t qstart = seg.qexon.start;
    uint32_t qend = seg.qexon.end;
    uint32_t len = qend - qstart;

    Cigar &cigar = match.align.cigar;

    if (k == 0 || k == (n - 1)) {
      cigar.add_operation(len, BAM_CSOFT_CLIP);
      match.prev_op = BAM_CSOFT_CLIP;
    } else {
      cigar.add_operation(len, BAM_CINS);
      match.prev_op = BAM_CINS;
    }
    match.total_operations += len;
    match.total_coverage += len;
    if (config.print) printf("adding ins len %d\n", len);
  }

  void ReadEvaluator::build_cigar_gap(Segment &seg, ExonChainMatch &match, 
                                      ReadEvaluationConfig &config) {
    uint32_t gstart = seg.gexon.start;
    uint32_t gend = seg.gexon.end;
    uint32_t len = gend - gstart;

    Cigar &cigar = match.align.cigar;

    cigar.add_operation(len, BAM_CDEL);
    match.prev_op = BAM_CDEL;
    match.total_operations += len;
    match.total_coverage += len;
    match.ref_consumed += len;
    if (config.print) printf("adding gap len %d\n", len);
  }

  void ReadEvaluator::build_cigar_clip(Segment &seg, ExonChainMatch &match, 
                                      ReadEvaluationConfig &config) {

    Cigar &cigar = match.align.cigar;
    for (int i = 0; i < seg.cigar.n_cigar; i++) {
      auto cig = seg.cigar.cigar[i];
      uint8_t op = cig & BAM_CIGAR_MASK;
      uint32_t len = cig >> BAM_CIGAR_SHIFT;
      cigar.add_operation(len, op);
      
      if (op == BAM_CMATCH_OVERRIDE 
        || op == BAM_CDEL_OVERRIDE) {
        match.ref_consumed += len;
      }
    }
    if (config.print) printf("Added clip cigar to match cigar\n");
    match.align.clip_score += seg.score;
  }

  void 
  ReadEvaluator::filter_by_similarity(std::vector<ExonChainMatch> &matches, 
                                      std::shared_ptr<g2tTree> g2t, ReadEvaluationConfig config) {
    if (LONG_READS) {
      for (auto it = matches.begin(); it != matches.end(); ) {
        auto &match = (*it);

        double similarity = (match.total_operations > 0)
          ? (match.total_coverage / match.total_operations)
          : 0.0;

        if (similarity > config.similarity_threshold) {
          double x = ((similarity - config.similarity_threshold) 
            / (1.0 - config.similarity_threshold));
          match.align.similarity_score = (x * x * static_cast<double>(match.junc_hits + 1)); // static_cast<double>(match.junc_hits + 1)
          
          if (config.print) {
            std::string tid_string = g2t->getTidName(match.tid);
            printf("x**3 = %f\n", x * x * x);
            printf("match.align.junc_hits + 1 = %f\n", static_cast<double>(match.junc_hits + 1));
            printf("tid = %s, similarity_score = %f\n", tid_string.c_str(), match.align.similarity_score);
          }

          ++it;
        } else {
          it = matches.erase(it);
        }
      }

      if (matches.empty()) return;
    }

    if (BRAMBLE_DEBUG) {
      for (auto it = matches.begin(); it != matches.end(); ) {
        auto &match = (*it);
        pos_t pos = match.align.strand == '+' ? 
          match.align.fwpos : match.align.rcpos;
        if (pos + match.ref_consumed > match.transcript_len) {
          std::string tid_string = g2t->getTidName(match.tid);
          printf("Error: read %s for transcript %s extends beyond transcript end\n",
            config.name.c_str(), tid_string.c_str());
          printf("transcript length = %d, reference consumed = %d, pos = %d\n",
            match.transcript_len, match.ref_consumed, pos);
          it = matches.erase(it);
        } else {
          ++it;
        }
      }

      if (matches.empty()) return;
    }
  }

  std::vector<ExonChainMatch> 
  ReadEvaluator::evaluate_exon_chains(CReadAln &read, 
                                      std::shared_ptr<g2tTree> g2t, 
                                      ReadEvaluationConfig config,
                                      uint8_t *seq, int seq_len) {

    uint32_t exon_count = read.segs.Count();
    refid_t refid = read.refid;
    config.name = read.brec->name();

    thread_local std::vector<ExonChainMatch> matches_by_strand;
    matches_by_strand.clear();

    bool has_left_clip = false;
    bool has_right_clip = false;
    uint32_t n_left_clip = 0;
    uint32_t n_right_clip = 0;

    bool failure = false;

    // if (config.name == "CHS.36908.12_ONT_simulated_read_266101") {
    //   config.print = true;
    // } else {
    //   config.print = false;
    // }
    // config.print = true;
    if (config.print) printf("read name = %s\n", config.name.c_str());

    thread_local unordered_map<tid_t, TidData> data;
    data.clear();

    if (LONG_READS) {
      get_clips(read, config, failure, has_left_clip, has_right_clip,
        n_left_clip, n_right_clip);
    }

    auto strands_to_check = get_strands_to_check(read);
    for (char strand : strands_to_check) {
      data.clear();
      failure = false;

      for (uint32_t j = 0; j < exon_count; j++) {
        if (config.print) printf("#exon is %d\n", j + 1);
        get_intervals(data, read, j, exon_count, config, g2t,
          refid, strand, has_left_clip, has_right_clip, failure);

        //if (!config.print) failure = true;
        if (failure) break;

      } // end of exon loop

      // try next strand
      if (failure) continue;

      if (LONG_READS && USE_FASTA) {
        for (auto &pair : data) {
          tid_t tid = pair.first;
          auto &td = pair.second;
          if (td.elim) continue;

          // Do soft clip rescue for long reads
          if (td.has_left_clip) {
            if (n_left_clip >= 5) {
              left_clip_rescue(td, strand, g2t, refid, 
                tid, n_left_clip, config, read, 
                seq, seq_len);
            } else {
              td.has_left_clip = false;
            }
          }

          if (td.has_right_clip) {
            if (n_right_clip >= 5) {
              right_clip_rescue(td, strand, g2t, refid,
                tid, n_right_clip, config, read, seq, seq_len);
            } else {
              td.has_right_clip = false;
            }
          }       
        }
      }

      // Calculate matches using segments
      // segments should now reflect ignored query exons, ignored guide exons,
      // and rectified soft clips on either side
      for (auto &pair : data) {
        tid_t tid = pair.first;
        auto &td = pair.second;

        if (td.elim) continue;

        if (config.print) {
          printf("\n=== SEGMENTS FOR TID %d (%s) ===\n", tid, g2t->getTidName(tid).c_str());
          std::vector<Segment> segs = td.segments;
          for (uint32_t k = 0; k < segs.size(); k++) {
              Segment *s = &segs[k];
              const char* status_str = "UNKNOWN";
              switch (s->status) {
                  case FIRST_EXON:   status_str = "FIRST";   break;
                  case MIDDLE_EXON:  status_str = "MIDDLE";  break;
                  case LAST_EXON:    status_str = "LAST";    break;
                  case ONLY_EXON:    status_str = "ONLY";    break;
                  case INS_EXON:     status_str = "INS";     break;
                  case GAP_EXON:     status_str = "GAP";     break;
                  case LEFTC_EXON:   status_str = "LEFTC";   break;
                  case RIGHTC_EXON:  status_str = "RIGHTC";  break;
              }
              printf("  [%u] status=%-8s has_qexon=%d has_gexon=%d is_small=%d",
                  k, status_str, s->has_qexon, s->has_gexon, s->is_small_exon);
              if (s->has_qexon)
                  printf(" qexon=[%u-%u]", s->qexon.start, s->qexon.end);
              if (s->has_gexon)
                  printf(" gexon=[%u-%u] exon_id=%u pos=%d", 
                      s->gexon.start, s->gexon.end, s->gexon.exon_id, s->gexon.pos);
              printf("\n");
          }
          printf("\n");
        }

        uint32_t n_segments = td.segments.size();
        bool match_created = false;

        uint32_t prev_s = 0;
        uint32_t prev_e = 0;

        uint32_t first_match_idx = -1;
        uint32_t last_match_idx = -1;

        for (uint32_t k = 0; k < n_segments; k++) {
          Segment &seg = td.segments[k];

          // if two read exons somehow mapped to same guide exon
          // just quit
          if (seg.has_gexon) {
            if (seg.gexon.start == prev_s 
              && seg.gexon.end == prev_e) {
              td.elim = true;
              break;
            }
          }

          prev_s = seg.gexon.start;
          prev_e = seg.gexon.end;

          // Create exon chain matches 
          if (!match_created) {
            create_match(td, seg.gexon, tid, strand);
            match_created = true;
            first_match_idx++;
            last_match_idx++;
            if (config.print) printf("match created, fwpos = %d, rcpos = %d\n", 
              td.match.align.fwpos, td.match.align.rcpos);
          }

          else if (match_created && seg.status != INS_EXON) {
            // update rcpos
            last_match_idx++;
            if (strand == '-') {
              td.match.align.rcpos = seg.gexon.pos;
            }
          }
        }

        bool first_match = false;
        bool last_match = false;

        // second pass — build CIGAR
        for (uint32_t k = 0; k < n_segments; k++) {
          Segment &seg = td.segments[k];
          if (td.elim) break;

          first_match = (k == first_match_idx);
          last_match = (k == last_match_idx);

          // build CIGAR as we process each exon
          bool is_match = (seg.status == FIRST_EXON || seg.status == MIDDLE_EXON
            || seg.status == LAST_EXON || seg.status == ONLY_EXON);
          bool is_ins = (seg.status == INS_EXON);
          bool is_gap = (seg.status == GAP_EXON);
          bool is_clip = (seg.status == LEFTC_EXON || seg.status == RIGHTC_EXON);

          if (is_match) {
            build_cigar_match(seg, td, td.match, first_match, 
              last_match, config);
            
          } else if (is_ins) {
            build_cigar_ins(seg, k, n_segments, td.match, config);
            td.match.junc_hits -= (k == 0 || k == n_segments - 1) ? 1 : 2;

          } else if (is_gap) {
            build_cigar_gap(seg, td.match, config);
            td.match.junc_hits -= 2;

          } else if (is_clip) {
            build_cigar_clip(seg, td.match, config);
          }
          
        } // end of segment loop

        // shouldn't happen
        if (td.match.junc_hits < 0) {
          td.match.junc_hits = 0;
        }

        // record match for this strand
        if (!td.elim) matches_by_strand.emplace_back(td.match);
      }

      if (config.print) {
        for (auto &match : matches_by_strand) {
          std::string tid_string = g2t->getTidName(match.tid);
          printf("tid = %d, %s\n", match.tid, tid_string.c_str());
          printf("fwpos = %d\n", match.align.fwpos);
          printf("rcpos = %d\n", match.align.rcpos);
          printf("IDEAL CIGAR: ");
          for (uint32_t i = 0; i < match.align.cigar.n_cigar; i++) {
            auto cig = match.align.cigar.cigar[i];
            printf("%u%c ", cig >> BAM_CIGAR_SHIFT, 
              "MIDNSHP=XB,./;"[cig & BAM_CIGAR_MASK]);
          }
          printf("\n\n");
        }
      }

    } // end of strands to check

    if (!matches_by_strand.empty()) {
      filter_by_similarity(matches_by_strand, g2t, config);
    }
    return matches_by_strand;
  }

  std::vector<ExonChainMatch>
  ShortReadEvaluator::evaluate(CReadAln &read,
                              std::shared_ptr<g2tTree> g2t,
                              uint8_t *seq,
                              int seq_len) {

    uint32_t max_clip = 2;
    uint32_t max_ins = 0;
    uint32_t max_gap = 0;
    uint32_t max_junc_gap = 0;
    float similarity_threshold = 0.95;

    ReadEvaluationConfig config = {
      max_clip,                       // max clip size
      max_ins,                        // max insertion to intervals
        // for when there exist no guides to explain a portion of an exon
      max_gap,                        // max gap in reference to intervals
      false,                          // ignore small exons?
      0,                              // small exon size
      max_junc_gap,                   // max junction gap
      similarity_threshold,
      false,                           // print debug statements
      ""
    }; 

    auto matches = evaluate_exon_chains(read, g2t, config, seq, seq_len);
    return matches;
  }

  std::vector<ExonChainMatch> 
  LongReadEvaluator::evaluate(CReadAln &read, 
                              std::shared_ptr<g2tTree> g2t,
                              uint8_t *seq,
                              int seq_len) {

    uint32_t max_clip = 40;
    uint32_t max_ins = 40;
    uint32_t max_gap = 40;
    uint32_t max_junc_gap = 40;
    float similarity_threshold = 0.60;
    
    ReadEvaluationConfig config = {
      max_clip,                       // max clip size
      max_ins,                        // max insertion to intervals
        // for when there exist no guides to explain a portion of an exon
      max_gap,                        // max gap in reference to intervals
      true,                           // ignore small exons?
      35,                             // small exon size
      max_junc_gap,                   // max junction gap
      similarity_threshold,
      false,                          // print debug statements
      ""
    };   

    auto matches = evaluate_exon_chains(read, g2t, config, seq, seq_len); 
    return matches;
  }

} // namespace bramble