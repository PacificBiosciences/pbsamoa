#ifndef PBSAMOA_IO_BAMMERGE_HPP
#define PBSAMOA_IO_BAMMERGE_HPP

#include <pbsamoa/io/BamSort.hpp>  // SortOrder

#include <array>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {

/// \brief Configuration for merging BAM files.
struct MergeConfig
{
    /// Sort order the inputs are sorted by (and the order of the merged output).
    /// std::nullopt auto-detects from the inputs' @HD SO: all-coordinate or
    /// all-queryname performs that sorted merge, all-unsorted concatenates, and a
    /// mix is rejected. An explicit value forces (and validates) that sorted order.
    std::optional<SortOrder> Order{};

    /// Two-char tag the inputs are sorted by; required iff Order == TAG.
    std::array<char, 2> Tag{};

    /// 0 -> auto (min(hardware_concurrency, 8)). Default size of both CPU worker
    /// pools on the sorted-merge path (input decompression and output compression);
    /// DecodeThreads / CompressThreads override each independently. Does not bound
    /// the per-input producer threads (one lightweight, I/O-bound thread per input).
    std::size_t NumThreads{0};

    /// Input-decompression pool size on the sorted-merge path. 0 -> inherit the
    /// resolved NumThreads. The dominant CPU consumer when inputs are highly
    /// compressed; raise it independently of the output compressor.
    std::size_t DecodeThreads{0};

    /// Output BGZF compression pool size on the sorted-merge path. 0 -> inherit the
    /// resolved NumThreads. Lower it (vs DecodeThreads) when the output level is
    /// cheap relative to input decompression.
    std::size_t CompressThreads{0};

    /// Output BGZF compression level [1, 12].
    int CompressionLevel{6};

    /// Total in-flight decompressed read-ahead budget across all inputs on the
    /// sorted-merge path (concat ignores it). A single cap, independent of input
    /// count; each source still keeps its head record available so a budget below
    /// one record cannot stall the merge. Default 768 MiB (matches sort's --memory).
    /// Peak resident RAM also includes the per-input current batch the consumer
    /// holds outside this budget (~= NumInputs x BatchBytes) plus the output
    /// writer's queue (WriterQueueCapacity blocks); size those knobs to bound them.
    ByteLimit ReadAheadMemory{std::size_t{768} * 1024 * 1024};

    /// Target framed-record payload per per-input handoff batch on the sorted-merge
    /// path. 0 -> library default (256 KiB). Each input holds one such batch outside
    /// ReadAheadMemory, so this sets the per-input working set above the budget;
    /// smaller trims that floor at the cost of more producer/consumer handoffs.
    std::size_t BatchBytes{0};

    /// Output BGZF writer input-queue depth (compressed-block slots) on the
    /// sorted-merge path. 0 -> library default (256). Bounds the writer's in-flight
    /// output memory; must be >= 1.
    std::size_t WriterQueueCapacity{0};

    /// Command line recorded in the appended @PG CL field (CLI passes argv).
    std::optional<std::string> CommandLine{};

    /// Force byte-level concatenation (BGZF block passthrough), ignoring input sort
    /// order. The output is marked SO:unsorted.
    bool Concat{false};

    /// When set, write a BAI index to this path, built on the fly during the merge.
    /// Requires a coordinate-sorted merge that takes the heap-merge path: BAI cannot
    /// index unsorted (--concat) or queryname/tag output, and the disjoint-chain
    /// passthrough emits verbatim BGZF blocks without decoding records, so any of those
    /// combined with BaiOutput is rejected (run bai-build on the output instead).
    std::optional<std::filesystem::path> BaiOutput{};
};

/// \brief Runtime profile of a merge for bottleneck diagnosis.
///
/// The k-way merge funnels every record through one consumer thread: Next()
/// (input) -> heap compare (cpu) -> writer.Write() (output), while producers
/// throttle on the read-ahead budget (memory). Each limiter is timed at its real
/// block point, sampled only on the slow path (an actual stall), so a smooth run
/// reports near-zero waits at no hot-path cost. The consumer's own compute is
/// measured directly (MergeCpuSeconds, a per-thread CPU clock) rather than inferred
/// as wall minus waits: the difference isolates CPU-starvation (the compress pool
/// oversubscribing cores) instead of mislabeling it as merge work. The CLI turns
/// these into a verdict.
///
/// Seconds fields are wall time. The heap-only fields are zero unless HeapMerge;
/// the concat / disjoint-passthrough paths are BGZF block copies with no decode or
/// compress pipeline, so only Wall / BytesIn / BytesOut are meaningful there.
struct MergeRuntimeStats
{
    bool HeapMerge{false};     ///< true iff the parallel k-way heap path ran
    double WallSeconds{0.0};   ///< wall time of the merge body
    std::int64_t BytesIn{0};   ///< summed compressed on-disk input bytes
    std::int64_t BytesOut{0};  ///< compressed on-disk output bytes

    double InputWaitSeconds{0.0};        ///< consumer starved for decoded records in Next()
    double OutputWaitSeconds{0.0};       ///< consumer spun on a full writer input queue
    double BudgetWaitSeconds{0.0};       ///< producers blocked on the read-ahead budget
    double InputIoReadSeconds{0.0};      ///< aggregate producer file-read time
    double InputDecompressSeconds{0.0};  ///< aggregate producer decompress time
    double OutputCompressSeconds{0.0};   ///< aggregate compress-pool time
    double OutputWriteSeconds{0.0};      ///< output file-write time
    double MergeCpuSeconds{0.0};         ///< measured CPU of the consumer (heap-merge) thread
    std::int64_t WriterStalls{0};        ///< IO writer found no ready block (compress-limited)
    std::int64_t PeakInFlightBytes{0};   ///< read-ahead high-water bytes
    std::int64_t BudgetBytes{0};         ///< read-ahead budget cap
    std::size_t DecodeThreads{0};        ///< resolved input-decompression pool size
    std::size_t CompressThreads{0};      ///< resolved output-compression pool size
};

/// \brief Outcome statistics from a merge.
struct MergeStats
{
    std::int64_t NumRecords{0};   ///< total records written
    std::size_t NumInputs{0};     ///< number of input files merged
    bool Passthrough{false};      ///< true iff a sorted merge emitted via verbatim BGZF
                                  ///< block copy (disjoint coordinate inputs)
    MergeRuntimeStats Runtime{};  ///< runtime profile for the CLI's bottleneck report
};

/// \brief Merge several BAM files into one.
///
/// All inputs must share an identical @SQ reference list; the merged header unions
/// their @RG / @PG / @CO records. Output is written atomically. The mode follows
/// \p config (see MergeConfig::Order / MergeConfig::Concat):
///   - Sorted merge: streaming k-way merge of inputs already sorted by the resolved
///     order. Equal-key records are emitted in input-file order (first input wins
///     ties), which is deterministic and independent of the thread count. When the
///     inputs of a coordinate merge form a disjoint chain (each file's coordinates
///     strictly precede the next's), the merge instead emits via verbatim BGZF block
///     passthrough — no recompression — and sets MergeStats::Passthrough.
///   - Concat: byte-level concatenation via BGZF block passthrough (no
///     decompress/recompress). NumRecords is reported as -1 (not counted).
///
/// \throws std::runtime_error on I/O errors, a missing input, mismatched/mixed sort
///         order, or incompatible reference lists (mirrors BamWriter / SortBam).
MergeStats MergeBam(const std::vector<std::filesystem::path>& inputs,
                    const std::filesystem::path& output, const MergeConfig& config = {});

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_IO_BAMMERGE_HPP
