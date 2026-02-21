#ifndef PBSAMOA_CORE_METRICS_HPP
#define PBSAMOA_CORE_METRICS_HPP

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {

/// \brief Snapshot of BGZF pipeline metrics (IO + decompression + output queue).
///
/// All fields are plain values captured at a point in time via relaxed atomics.
/// Diff two snapshots to compute per-second rates.
struct BgzfMetrics
{
    // --- IO stage throughput ---
    std::uint64_t BytesRead{0};          ///< compressed bytes read from disk
    std::uint64_t BlocksRead{0};         ///< BGZF blocks submitted to pool
    std::uint64_t BytesDecompressed{0};  ///< decompressed bytes produced

    // --- Decompression pool (forwarded from ThreadPool::Metrics) ---
    std::size_t PoolQueueDepth{0};         ///< current work queue depth
    std::size_t PoolPeakQueueDepth{0};     ///< max work queue depth observed
    std::size_t PoolActiveWorkers{0};      ///< currently executing decompression tasks
    std::size_t PoolPeakActiveWorkers{0};  ///< max concurrent active workers observed
    std::size_t PoolResultQueueDepth{0};   ///< pending results for consumer thread
    std::size_t PoolPeakResultQueueDepth{0};

    // --- Pipeline output SPSC queue ---
    std::uint64_t RecordsProduced{0};  ///< records pushed by consumer thread
    std::uint64_t RecordsConsumed{0};  ///< records popped by caller
    /// Queue depth = RecordsProduced - RecordsConsumed

    // --- Stall counters (direct bottleneck indicators) ---
    std::uint64_t IoStalls{0};        ///< IO thread found pool queue at capacity before Submit
    std::uint64_t ConsumerStalls{0};  ///< consumer thread failed SPSC try_push (queue full)
    std::uint64_t ReaderStalls{0};    ///< caller entered cv.wait slow path (queue empty)

    // --- Cumulative stage timing (nanoseconds, per-batch granularity) ---
    std::uint64_t IoReadNs{0};       ///< time in file reads
    std::uint64_t DecompressNs{0};   ///< time in libdeflate (summed across all workers)
    std::uint64_t RecordParseNs{0};  ///< time parsing records in consumer thread
};

/// \brief Snapshot of BAM record decode metrics (BamRecordReader layer).
struct DecodeMetrics
{
    // --- Decode pool (forwarded from ThreadPool::Metrics) ---
    std::size_t PoolQueueDepth{0};
    std::size_t PoolPeakQueueDepth{0};
    std::size_t PoolActiveWorkers{0};
    std::size_t PoolPeakActiveWorkers{0};

    // --- Throughput ---
    std::uint64_t BatchesDecoded{0};
    std::uint64_t RecordsDecoded{0};

    // --- Output SPSC queue ---
    std::uint64_t RecordsProduced{0};
    std::uint64_t RecordsConsumed{0};
    /// Queue depth = RecordsProduced - RecordsConsumed

    // --- Stall counters ---
    std::uint64_t ProducerStalls{0};  ///< producer blocked by full output queue
    std::uint64_t ConsumerStalls{0};  ///< caller yield-spun on empty output queue

    // --- Timing (nanoseconds) ---
    std::uint64_t DecodeNs{0};     ///< cumulative decode time (summed across all workers)
    std::uint64_t BatchReadNs{0};  ///< time spent in ReadBatch from underlying reader
};

/// \brief Composite snapshot returned by BamRecordReader::GetMetrics().
struct ReaderMetrics
{
    BgzfMetrics Bgzf;
    DecodeMetrics Decode;
    std::uint64_t TotalRecordsRead{0};  ///< total records returned to caller
    bool ParallelBgzf{false};           ///< true if BgzfWorkers > 0
    bool ParallelDecode{false};         ///< true if DecodeWorkers > 0
};

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_CORE_METRICS_HPP
