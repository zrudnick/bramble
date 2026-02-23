#include <numeric>
#include <algorithm>
#include <iostream>
#include <memory>
#include <vector>
#include <atomic>

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
#include "bramble.h"
#include "evaluate.h"

extern bool VERBOSE;
extern bool BRAMBLE_DEBUG;

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
extern std::atomic<uint8_t> threads_waiting; // Idle worker threads
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
      are_threads_waiting = 
        (threads_waiting.load(std::memory_order_relaxed) > 0);

      if (are_threads_waiting) {
        have_bundles.notify_all();
        current_thread::sleep_for(1);

        are_threads_waiting = 
          (threads_waiting.load(std::memory_order_relaxed) > 0);

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
    int threads = 
      threads_waiting.load(std::memory_order_relaxed);
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
      threads_waiting.fetch_add(1, std::memory_order_relaxed);
      queue_mutex.unlock();
      wait_mutex.lock();
      have_threads.notify_one(); // in case main thread is waiting
      wait_mutex.unlock();
      current_thread::yield();
      queue_mutex.lock();
      while (bundle_work && bundle_queue->Count() == 0) {
        // Unlocks queue_mutex and wait until notified
        // When notified, locks queue_mutex and resume
        have_bundles.wait(queue_mutex);
      }

      if (threads_waiting.load(std::memory_order_relaxed) > 0)
        threads_waiting.fetch_sub(1, std::memory_order_relaxed);

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

  // Process reads in previous bundle
  void push_bundle(BundleData *bundle, GPVec<BundleData> *bundle_queue) {

    if (bundle->reads.size() > 0) {
      bundle->getReady();

#ifndef NOTHREADS

      // Push this in the bundle queue where it'll be picked up by the threads
      int queue_count = 0;

      queue_mutex.lock();
      bundle_queue->Push(bundle);
      bundle_work |= 0x02; // set bit 1 to 1
      queue_count = bundle_queue->Count();
      queue_mutex.unlock();

      wait_mutex.lock();
      while (threads_waiting.load(std::memory_order_relaxed) == 0) {
        have_threads.wait(wait_mutex);
      }
      wait_mutex.unlock();
      have_bundles.notify_one();

      current_thread::yield();

      queue_mutex.lock();
      while (bundle_queue->Count() == queue_count) {
        queue_mutex.unlock();
        have_bundles.notify_one();
        current_thread::yield();
        queue_mutex.lock();
      }
      queue_mutex.unlock();

#else // NOTHREADS = true

      // Just process single bundle
      process_bundle(bundle, io);
#endif
    }

    // Clear bundle (no more alignments)
    else {
      #ifndef NOTHREADS
        data_mutex.lock();
      #endif

      bundle->Clear();

      #ifndef NOTHREADS
        clear_data_pool.Push(bundle->idx);
        data_mutex.unlock();
      #endif
    }
  }
}