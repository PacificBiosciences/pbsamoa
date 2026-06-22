#ifndef PBSAMOA_IO_MERGEREADAHEAD_HPP
#define PBSAMOA_IO_MERGEREADAHEAD_HPP

#include <pbsamoa/io/BamRawReader.hpp>  // ByteLimit

#include <filesystem>
#include <memory>
#include <optional>
#include <span>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {

/// \brief Bottleneck-diagnosis counters sampled across a read-ahead's lifetime.
///
/// All times are wall nanoseconds accumulated only on the slow path (an actual
/// block), so a never-starved pipeline reports zeros at no hot-path cost. Read
/// after the merge drains (all producers finished); fields are relaxed atomics.
struct MergeReadAheadStats
{
    /// Wall time the consumer spent blocked in Next() waiting for a batch. High
    /// => decode/read cannot feed the merge fast enough (input-bound).
    std::uint64_t ConsumerInputWaitNs{0};

    /// Wall time producers spent blocked on a full read-ahead budget (summed
    /// across producers). High with PeakInFlightBytes ~= BudgetBytes => the
    /// budget caps throughput (memory-bound; raise --memory).
    std::uint64_t ProducerBudgetWaitNs{0};

    /// Aggregate producer time in file reads vs parallel decompress (summed
    /// across producers). Their ratio splits an input-bound run into disk-read
    /// vs decode-CPU limited.
    std::uint64_t IoReadNs{0};
    std::uint64_t DecompressNs{0};

    std::size_t PeakInFlightBytes{0};  ///< high-water summed in-flight batch bytes
    std::size_t BudgetBytes{0};        ///< configured read-ahead budget
};

/// \brief Parallel, memory-bounded read-ahead over several BAM inputs for the
/// sorted k-way merge.
///
/// Each input gets a background producer that reads raw BGZF blocks and
/// decompresses them on a shared thread pool, framing records into per-source
/// batches handed to the consumer as views into an internal buffer (no per-record
/// copy or allocation). A single \p budget caps the summed in-flight decompressed
/// bytes across all sources; producers sleep on a condition variable when the
/// budget (or a queue) is full rather than busy-spinning. Each source always keeps
/// at least its head batch available, so the merge cannot stall even under a budget
/// smaller than one batch. Records are delivered per source strictly in file order,
/// so the merge heap's input-order tie-breaking — and the merged output — are
/// independent of thread count, decode timing, and budget.
class MergeReadAhead
{
public:
    /// \param[in] inputs        BAM files; records are read after each file's header.
    /// \param[in] decodeWorkers shared decompression pool size (clamped to >= 1).
    /// \param[in] budget        total in-flight decompressed read-ahead bytes.
    /// \param[in] batchBytes    target framed-record payload per per-source handoff
    ///                          batch; 0 -> default (256 KiB). Each source holds one
    ///                          current batch outside \p budget, so this sets the
    ///                          per-source working set floor above the budget.
    /// \throws std::runtime_error if an input cannot be opened.
    MergeReadAhead(std::span<const std::filesystem::path> inputs, std::size_t decodeWorkers,
                   ByteLimit budget, std::size_t batchBytes = 0);
    ~MergeReadAhead();

    MergeReadAhead(const MergeReadAhead&) = delete;
    MergeReadAhead& operator=(const MergeReadAhead&) = delete;
    MergeReadAhead(MergeReadAhead&&) = delete;
    MergeReadAhead& operator=(MergeReadAhead&&) = delete;

    std::size_t NumSources() const noexcept;

    /// \brief Next record from source \p sourceIndex in file order, as a view into
    /// an internal buffer (the BAM record payload after block_size, starting at
    /// refID — the same span RawRecord wraps).
    /// \returns the record bytes, or std::nullopt at that source's end.
    /// \warning The span is valid until the next Next(sourceIndex) call; copy out
    ///          anything you need to retain past that. It outlives a single
    ///          BamWriter::Write(span), which copies synchronously.
    /// \throws the producer's captured exception on a read/decompress error.
    [[nodiscard]] std::optional<std::span<const std::byte>> Next(std::size_t sourceIndex);

    /// \brief Bottleneck counters for the merge runtime report. Call after the
    /// merge has drained every source (the consumer saw EOF on all of them).
    [[nodiscard]] MergeReadAheadStats Stats() const noexcept;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_IO_MERGEREADAHEAD_HPP
