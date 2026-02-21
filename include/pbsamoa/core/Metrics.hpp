#ifndef PBSAMOA_CORE_METRICS_HPP
#define PBSAMOA_CORE_METRICS_HPP

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {

/// \brief Snapshot of ThreadPool metrics (forwarded from ThreadPool::Metrics).
struct PoolMetrics
{
    std::size_t QueueDepth{0};         ///< current work queue depth
    std::size_t PeakQueueDepth{0};     ///< max work queue depth observed
    std::size_t ActiveWorkers{0};      ///< currently executing tasks
    std::size_t PeakActiveWorkers{0};  ///< max concurrent active workers observed
    std::size_t ResultQueueDepth{0};   ///< pending results for consumer thread
    std::size_t PeakResultQueueDepth{0};
};

/// \brief Snapshot of BGZF pipeline metrics (IO + decompression + output
/// queue).
///
/// All fields are plain values captured at a point in time via relaxed atomics.
/// Diff two snapshots to compute per-second rates.
struct BgzfMetrics
{
    // --- IO stage throughput ---
    std::uint64_t BytesRead{0};          ///< compressed bytes read from disk
    std::uint64_t BlocksRead{0};         ///< BGZF blocks submitted to pool
    std::uint64_t BytesDecompressed{0};  ///< decompressed bytes produced

    // --- Decompression pool ---
    PoolMetrics Pool;

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
    // --- Decode pool ---
    PoolMetrics Pool;

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

/// \brief Snapshot of BGZF writer pipeline metrics.
struct BgzfWriteMetrics
{
    // --- Stall counters ---
    std::uint64_t CallerStalls{0};  ///< caller blocked on full input queue
    std::uint64_t PackerStalls{0};  ///< packer found input queue empty
    std::uint64_t WriterStalls{0};  ///< IO writer found no ready results

    // --- Throughput ---
    std::uint64_t BytesCompressed{0};  ///< compressed payload bytes produced
    std::uint64_t BlocksWritten{0};    ///< BGZF blocks written

    // --- Timing (nanoseconds) ---
    std::uint64_t CompressNs{0};  ///< cumulative libdeflate time
    std::uint64_t IoWriteNs{0};   ///< cumulative file write time
    std::uint64_t CallbackNs{0};  ///< cumulative callback dispatch time

    // --- Compression pool ---
    PoolMetrics Pool;
};

/// \brief Composite writer metrics returned by BamWriter::GetMetrics().
struct WriterMetrics
{
    BgzfWriteMetrics Bgzf;
    std::uint64_t TotalRecordsWritten{0};
};

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_CORE_METRICS_HPP
