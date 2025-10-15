
#pragma once

namespace bramble {

  struct BundleData;
  struct BamIO;

  bool has_more_bundles();

  void wait_for_bundles();

  void process_bundle(BundleData *bundle, BamIO *io);

  bool no_threads_waiting();

  void worker_thread(GThreadData &td);

  int wait_for_data(BundleData *bundles);
}