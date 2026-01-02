#include <numeric>
#include <algorithm>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#ifndef NOTHREADS
#include "GThreads.h"
#endif

#include "GArgs.h"
#include "GBitVec.h"
#include "GSam.h"
#include "GStr.h"
#include "gff.h"
#include "time.h"

#include "types.h"
#include "bundles.h"
#include "bramble.h"
#include "evaluate.h"

extern bool VERBOSE;
extern bool DEBUG;

extern uint8_t n_threads;    // Threads, -p

#ifndef NOTHREADS
// Threading: single producer, multiple consumers
// Main thread/program is always loading the producer

extern GMutex data_mutex; // Manage availability of data records ready to be loaded by
                   // main thread
extern GVec<int> clear_data_pool; // Indices of data bundles cleared for loading by
                           // main thread (clear data pool)
extern GConditionVar
    have_bundles; // Will notify a thread that a bundle was loaded in the ready
                  // queue (or that no more bundles are coming)
extern char bundle_work; // Bit 0 set if bundles are still being prepared (BAM file not exhausted
       // yet) Bit 1 set if there are Bundles ready in the queue

extern GMutex wait_mutex;       // Controls threads_waiting (idle threads counter)
extern uint8_t threads_waiting; // Idle worker threads
extern GConditionVar
    have_threads; // Will notify the bundle loader when a thread
                  // Is available to process the currently loaded bundle

extern GConditionVar have_clear; // Will notify when bundle buf space available
extern GMutex queue_mutex;       // Controls bundle_queue and bundles access
extern GFastMutex log_mutex;     // Only when verbose - to avoid mangling the log output
extern GFastMutex reading_mutex;
extern GFastMutex bam_io_mutex;  // Protects BAM io
#endif

extern bool no_more_bundles;

namespace bramble {

  // Check if there are more bundles
  bool has_more_bundles() {
#ifndef NOTHREADS
    GLockGuard<GFastMutex> lock(reading_mutex);
#endif
    return !no_more_bundles;
  }

  // Wait for threads and set bundle status
  void wait_for_bundles() {

#ifndef NOTHREADS
    reading_mutex.lock();
    no_more_bundles = true;
    reading_mutex.unlock();

    queue_mutex.lock();
    bundle_work &= ~(int)0x01; // clear bit 0;
    queue_mutex.unlock();

    bool are_threads_waiting = true;
    while (are_threads_waiting) {
      wait_mutex.lock();
      are_threads_waiting = (threads_waiting > 0);
      wait_mutex.unlock();

      if (are_threads_waiting) {
        have_bundles.notify_all();
        current_thread::sleep_for(1);

        wait_mutex.lock();
        are_threads_waiting = (threads_waiting > 0);
        wait_mutex.unlock();

        current_thread::sleep_for(1);
      }
    }
#else
    no_more_bundles = true;
#endif
  }

  // Process current bundle
  void process_bundle(BundleData *bundle, BamIO *io) {
    convert_reads(bundle->reads, bundle->g2t, bundle->evaluator, io);
    bundle->Clear();
  }

  // Check that there aren't any threads waiting
  bool no_threads_waiting() {
    wait_mutex.lock();
    int threads = threads_waiting;
    wait_mutex.unlock();
    return (threads < 1);
  }

  // Worker thread waits for incoming bundle, calls process_bundle
  void worker_thread(GThreadData &td) {
    WorkerArgs *args = (WorkerArgs *)(td.udata);
    GPVec<BundleData> *bundle_queue = args->bundle_queue;
    BamIO *io = args->io;

    // Wait for a ready bundle in the queue
    queue_mutex.lock(); // enter wait-for-notification loop

    while (bundle_work) {
      wait_mutex.lock();
      threads_waiting++;
      queue_mutex.unlock();
      wait_mutex.unlock();
      have_threads.notify_one(); // in case main thread is waiting
      current_thread::yield();
      queue_mutex.lock();
      while (bundle_work && bundle_queue->Count() == 0) {
        // Unlocks queue_mutex and wait until notified
        // When notified, locks queue_mutex and resume
        have_bundles.wait(queue_mutex);
      }

      wait_mutex.lock();
      if (threads_waiting > 0)
      threads_waiting--;
      wait_mutex.unlock();

      BundleData *readyBundle = NULL;
      if ((bundle_work & 0x02) != 0) {
        readyBundle = bundle_queue->Pop();

        if (readyBundle != NULL) {
          if (bundle_queue->Count() == 0)
          bundle_work &= ~(int)0x02; // clear bit 1 (queue is empty)

          queue_mutex.unlock();
          process_bundle(readyBundle, io);

          data_mutex.lock();
          clear_data_pool.Push(readyBundle->idx);
          data_mutex.unlock();

          have_clear.notify_one(); // inform main thread
          current_thread::yield();

          queue_mutex.lock();
        }
      }
    }
    queue_mutex.unlock();
  }

  // Prepare the next available bundle slot for loading
  int wait_for_data(BundleData *bundles) {
    int idx = -1;

    data_mutex.lock();
    while (clear_data_pool.Count() == 0) {
      have_clear.wait(data_mutex);
    }
    idx = clear_data_pool.Pop();
    if (idx >= 0)
      bundles[idx].status = BundleStatus::BUNDLE_STATUS_LOADING;
    data_mutex.unlock();

    return idx;
  }
}