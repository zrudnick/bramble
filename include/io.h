
#pragma once

#include "GSam.h"
#include "GStr.h"

namespace bramble {

  struct BamIO {
  protected:
    std::unique_ptr<GSamReader> reader{nullptr};
    std::unique_ptr<GSamWriter> writer{nullptr};
    GSamRecord *current_record = nullptr;
    GStr input_bam;
    GStr output_bam;
    GStr header_sam;

  public:
    BamIO(const GStr &in_bam, const GStr &out_bam, const GStr &sam_header)
        : input_bam(in_bam), output_bam(out_bam), header_sam(sam_header) {}

    void start() {
      reader = std::make_unique<GSamReader>(input_bam.chars());
      writer = std::make_unique<GSamWriter>(output_bam.chars(), header_sam.chars());

      GSamRecord *next_rec = reader->next();
      if (next_rec) {
        current_record = next_rec;
      }
    }

    GSamRecord *next() {
      GSamRecord *result = current_record;

      current_record = reader->next(); // null at EOF

      return result;
    }

    void write(GSamRecord *rec) {
      if (writer && rec != nullptr) {
        writer->write(rec);
      }
    }

    int32_t get_tid(const char *transcript_id) {
      return writer->get_tid(transcript_id);
    }

    void stop() {
      if (current_record != nullptr) {
        delete current_record;
        current_record = nullptr;
      }
    }
  };

}