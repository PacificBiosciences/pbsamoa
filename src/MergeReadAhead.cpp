#include <pbsamoa/io/MergeReadAhead.hpp>

#include "BinaryUtils.hpp"  // ComputeHeaderSize, MAX_DECOMPRESSED_BLOCK_SIZE

#include <pbsamoa/core/Bgzf.hpp>

#include <pbcopper/parallel/ThreadPool.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <deque>
#include <exception>
#include <filesystem>
#include <format>
#include <fstream>
#include <memory>
#include <mutex>
#include <optional>
#include <span>
#include <stdexcept>
#include <stop_token>
#include <thread>
#include <utility>
#include <vector>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {

namespace {

using DecodePool = PacBio::Parallel::ThreadPool<>;

/// Compressed BGZF blocks read per IO step before a parallel decompress dispatch.
/// Matches BgzfReader's batch granularity so one large input saturates the pool.
constexpr std::int32_t MERGE_BLOCKS_PER_BATCH{32};

/// Target framed-record payload per handoff batch. Records are handed to the merge
/// consumer in batches of roughly this size so the per-record cost is an index
/// advance, not a mutex + condition-variable round trip. The consumer holds one
/// current batch per source outside the budget, so this also sets the minimum
/// working set above the budget.
constexpr std::size_t BATCH_TARGET_BYTES{std::size_t{256} * 1024};

/// Per-record bookkeeping slack added to the raw byte count for the budget meter.
constexpr std::size_t RECORD_SLACK{64};

std::uint64_t ElapsedNs(std::chrono::steady_clock::time_point start)
{
    return static_cast<std::uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(
                                          std::chrono::steady_clock::now() - start)
                                          .count());
}

}  // namespace

struct MergeReadAhead::Impl
{
    /// One framed record's location within a batch buffer (payload after the
    /// 4-byte block_size, i.e. the span RawRecord wraps).
    struct Extent
    {
        std::uint32_t Offset;
        std::uint32_t Length;
    };

    /// A run of framed records handed to the consumer as a unit. Owns the
    /// decompressed bytes the records view into (moved in, not copied per record).
    struct Batch
    {
        std::vector<std::byte> Buffer;
        std::vector<Extent> Extents;
        std::size_t Bytes{0};
    };

    struct PerSource
    {
        std::filesystem::path Path;

        std::mutex QueueMutex;
        std::condition_variable DataCv;  // consumer waits for a batch / EOF / error
        std::deque<Batch> Queue;
        std::atomic<std::size_t> QueuedBatches{0};  // read lock-free by the budget predicate
        bool Done{false};
        std::exception_ptr Error{};

        // Consumer-side cursor — touched only by the single consumer thread.
        Batch Current{};
        std::size_t CurrentPos{0};

        std::jthread Producer;
    };

    // Pool is declared before Sources so producers (which use the pool) are joined
    // on destruction before the pool is torn down.
    std::shared_ptr<DecodePool> Pool;

    std::mutex BudgetMutex;
    std::condition_variable SpaceCv;  // producers wait for budget room
    std::size_t InFlight{0};          // summed queued batch bytes across all sources
    std::size_t BudgetBytes{0};
    std::size_t BatchTargetBytes{BATCH_TARGET_BYTES};  // per-source flush threshold
    bool Stopping{false};  // set under BudgetMutex to wake producers on teardown

    // Bottleneck-diagnosis counters (see MergeReadAheadStats). All accumulate only
    // on the slow path, so a never-starved pipeline pays nothing. PeakInFlight is
    // updated under BudgetMutex; the time counters are relaxed across threads.
    std::atomic<std::uint64_t> ConsumerInputWaitNs{0};
    std::atomic<std::uint64_t> ProducerBudgetWaitNs{0};
    std::atomic<std::uint64_t> IoReadNs{0};
    std::atomic<std::uint64_t> DecompressNs{0};
    std::atomic<std::size_t> PeakInFlightBytes{0};

    std::vector<std::unique_ptr<PerSource>> Sources;

    Impl(std::span<const std::filesystem::path> inputs, std::size_t decodeWorkers, ByteLimit budget,
         std::size_t batchBytes)
        : Pool{std::make_shared<DecodePool>(DecodePool::Config{
              .NumThreads = std::max<std::size_t>(decodeWorkers, 1),
              .QueueMultiplier = 3,
          })}
        , BudgetBytes{budget.Value()}
        , BatchTargetBytes{(batchBytes != 0) ? batchBytes : BATCH_TARGET_BYTES}
    {
        Sources.reserve(std::size(inputs));
        for (const std::filesystem::path& input : inputs) {
            auto source{std::make_unique<PerSource>()};
            source->Path = input;
            // Probe-open so a missing/unreadable input fails in the constructor (the
            // same place the old synchronous path failed) rather than asynchronously.
            const std::ifstream probe{input, std::ios::binary};
            if (!probe) {
                throw std::runtime_error{
                    std::format("MergeReadAhead: cannot open input {}", input.string())};
            }
            Sources.push_back(std::move(source));
        }
        for (std::size_t i{0}; i < std::size(Sources); ++i) {
            Sources[i]->Producer =
                std::jthread{[this, i](std::stop_token stopToken) { ProducerLoop(stopToken, i); }};
        }
    }

    Impl(const Impl&) = delete;
    Impl& operator=(const Impl&) = delete;
    Impl(Impl&&) = delete;
    Impl& operator=(Impl&&) = delete;

    ~Impl()
    {
        {
            const std::lock_guard lock{BudgetMutex};
            Stopping = true;
        }
        SpaceCv.notify_all();
        for (const std::unique_ptr<PerSource>& source : Sources) {
            source->Producer.request_stop();
        }
        // Sources (jthreads) join on member destruction, before Pool is destroyed.
    }

    /// Decompress \p blocks in parallel on the shared pool and append the results,
    /// in block order, to \p carry. Throws on a decompression failure.
    void DecompressInto(const std::vector<std::vector<std::byte>>& blocks,
                        const std::vector<BgzfBlockInfo>& infos, std::vector<std::byte>& carry)
    {
        const std::int32_t n{static_cast<std::int32_t>(std::size(blocks))};
        std::vector<std::vector<std::byte>> outputs(std::size(blocks));
        const auto worker{[&](std::int32_t i) {
            std::vector<std::byte> buffer(MAX_DECOMPRESSED_BLOCK_SIZE);
            const std::optional<std::size_t> produced{DecompressBgzfBlock(
                blocks[static_cast<std::size_t>(i)], infos[static_cast<std::size_t>(i)], buffer)};
            if (!produced) {
                throw std::runtime_error{"MergeReadAhead: BGZF decompression failed"};
            }
            buffer.resize(*produced);
            outputs[static_cast<std::size_t>(i)] = std::move(buffer);
        }};
        PacBio::Parallel::Dispatch(Pool, worker, n);
        for (std::vector<std::byte>& output : outputs) {
            carry.insert(std::end(carry), std::begin(output), std::end(output));
        }
    }

    /// Block on the budget, then publish \p batch to \p source. Returns false if
    /// teardown was requested while waiting.
    bool PublishBatch(PerSource& source, Batch&& batch)
    {
        {
            std::unique_lock budgetLock{BudgetMutex};
            // Always admit a source's first queued batch (QueuedBatches == 0) so the
            // merge can make progress even under a budget smaller than one batch.
            const auto admit{[&] {
                return Stopping || (source.QueuedBatches.load() == 0) ||
                       ((InFlight + batch.Bytes) <= BudgetBytes);
            }};
            // Time only a genuine budget-full block (admit() already false and not a
            // teardown / first-batch admit) so it attributes cleanly to memory.
            if (!admit()) {
                const auto blockStart{std::chrono::steady_clock::now()};
                SpaceCv.wait(budgetLock, admit);
                ProducerBudgetWaitNs.fetch_add(ElapsedNs(blockStart), std::memory_order_relaxed);
            }
            if (Stopping) {
                return false;
            }
            InFlight += batch.Bytes;
            if (InFlight > PeakInFlightBytes.load(std::memory_order_relaxed)) {
                PeakInFlightBytes.store(InFlight, std::memory_order_relaxed);
            }
        }
        {
            const std::lock_guard queueLock{source.QueueMutex};
            source.Queue.push_back(std::move(batch));
            source.QueuedBatches.fetch_add(1);
        }
        source.DataCv.notify_one();
        return true;
    }

    void ProducerLoop(std::stop_token stopToken, std::size_t index)
    {
        PerSource& source{*Sources[index]};
        try {
            std::ifstream in{source.Path, std::ios::binary};
            if (!in) {
                throw std::runtime_error{
                    std::format("MergeReadAhead: cannot open input {}", source.Path.string())};
            }

            std::vector<std::byte> carry{};  // working buffer of decompressed bytes
            std::vector<Extent> extents{};   // framed records, offsets into carry
            std::size_t pendingBytes{0};
            std::size_t framedEnd{0};  // bytes of carry already framed
            bool headerStripped{false};

            // Move the framed prefix [0, framedEnd) of carry into a batch (no
            // per-record copy), keep the unframed tail, and publish. Returns false
            // on teardown.
            const auto flush{[&]() -> bool {
                std::vector<std::byte> tail(
                    std::begin(carry) + static_cast<std::ptrdiff_t>(framedEnd), std::end(carry));
                Batch batch{};
                batch.Buffer = std::move(carry);
                batch.Buffer.resize(framedEnd);
                batch.Extents = std::move(extents);
                batch.Bytes = pendingBytes;

                carry = std::move(tail);
                extents.clear();
                pendingBytes = 0;
                framedEnd = 0;
                return PublishBatch(source, std::move(batch));
            }};

            while (!stopToken.stop_requested()) {
                std::vector<std::vector<std::byte>> blocks{};
                std::vector<BgzfBlockInfo> infos{};
                blocks.reserve(MERGE_BLOCKS_PER_BATCH);
                infos.reserve(MERGE_BLOCKS_PER_BATCH);
                bool eof{false};

                const auto ioStart{std::chrono::steady_clock::now()};
                for (std::int32_t b{0}; b < MERGE_BLOCKS_PER_BATCH; ++b) {
                    std::array<std::byte, 18> blockHeader{};
                    in.read(reinterpret_cast<char*>(std::data(blockHeader)), 18);
                    const std::streamsize headerRead{in.gcount()};
                    if (headerRead == 0) {
                        eof = true;
                        break;
                    }
                    if (headerRead < 18) {
                        throw std::runtime_error{std::format(
                            "MergeReadAhead: truncated BGZF header in {}", source.Path.string())};
                    }
                    const std::optional<BgzfBlockInfo> info{ParseBgzfBlockHeader(blockHeader)};
                    if (!info) {
                        throw std::runtime_error{std::format(
                            "MergeReadAhead: malformed BGZF block in {}", source.Path.string())};
                    }
                    std::vector<std::byte> block(info->blockSize);
                    std::ranges::copy(blockHeader, std::begin(block));
                    const std::streamsize bodyBytes{static_cast<std::streamsize>(info->blockSize) -
                                                    18};
                    in.read(reinterpret_cast<char*>(std::data(block) + 18), bodyBytes);
                    if (in.gcount() < bodyBytes) {
                        throw std::runtime_error{std::format(
                            "MergeReadAhead: truncated BGZF block in {}", source.Path.string())};
                    }
                    if (IsBgzfEofMarker(block)) {
                        eof = true;
                        break;
                    }
                    blocks.push_back(std::move(block));
                    infos.push_back(*info);
                }
                IoReadNs.fetch_add(ElapsedNs(ioStart), std::memory_order_relaxed);

                if (!blocks.empty()) {
                    const auto decStart{std::chrono::steady_clock::now()};
                    DecompressInto(blocks, infos, carry);
                    DecompressNs.fetch_add(ElapsedNs(decStart), std::memory_order_relaxed);

                    if (!headerStripped) {
                        const std::size_t headerLen{ComputeHeaderSize(carry)};
                        if (headerLen == 0) {
                            continue;  // header spans into a later block
                        }
                        carry.erase(std::begin(carry),
                                    std::begin(carry) + static_cast<std::ptrdiff_t>(headerLen));
                        headerStripped = true;
                    }

                    // Frame complete records from where we left off into extents.
                    std::size_t pos{framedEnd};
                    while ((std::size(carry) - pos) >= 4) {
                        const std::uint32_t recordLen{ReadU32LE(std::data(carry) + pos)};
                        if ((std::size(carry) - pos - 4) < recordLen) {
                            break;  // record continues into a later block
                        }
                        extents.push_back(Extent{.Offset = static_cast<std::uint32_t>(pos + 4),
                                                 .Length = recordLen});
                        pendingBytes += recordLen + RECORD_SLACK;
                        pos += 4 + recordLen;
                    }
                    framedEnd = pos;

                    if (pendingBytes >= BatchTargetBytes) {
                        if (!flush()) {
                            return;  // teardown
                        }
                    }
                }

                if (eof) {
                    break;
                }
            }

            if (!extents.empty()) {
                if (!flush()) {
                    return;  // teardown
                }
            }
            if (headerStripped && !carry.empty()) {
                throw std::runtime_error{std::format(
                    "MergeReadAhead: truncated record region in {}", source.Path.string())};
            }

            {
                const std::lock_guard lock{source.QueueMutex};
                source.Done = true;
            }
            source.DataCv.notify_all();
        } catch (...) {
            {
                const std::lock_guard lock{source.QueueMutex};
                source.Error = std::current_exception();
                source.Done = true;
            }
            source.DataCv.notify_all();
        }
    }

    std::optional<std::span<const std::byte>> Next(std::size_t index)
    {
        PerSource& source{*Sources[index]};

        // Fast path: hand out the next record view from the current batch lock-free.
        if (source.CurrentPos >= std::size(source.Current.Extents)) {
            // Current batch drained — fetch the next one.
            Batch next{};
            {
                std::unique_lock queueLock{source.QueueMutex};
                const auto ready{
                    [&] { return !source.Queue.empty() || source.Done || source.Error; }};
                // Time only a genuine starvation block, so it attributes to input.
                if (!ready()) {
                    const auto blockStart{std::chrono::steady_clock::now()};
                    source.DataCv.wait(queueLock, ready);
                    ConsumerInputWaitNs.fetch_add(ElapsedNs(blockStart), std::memory_order_relaxed);
                }
                if (source.Queue.empty()) {
                    if (source.Error) {
                        const std::exception_ptr error{source.Error};
                        queueLock.unlock();
                        std::rethrow_exception(error);
                    }
                    return std::nullopt;  // clean EOF
                }
                next = std::move(source.Queue.front());
                source.Queue.pop_front();
            }

            source.QueuedBatches.fetch_sub(1);
            {
                const std::lock_guard budgetLock{BudgetMutex};
                InFlight -= next.Bytes;
            }
            SpaceCv.notify_all();

            source.Current = std::move(next);
            source.CurrentPos = 0;
        }

        const Extent extent{source.Current.Extents[source.CurrentPos++]};
        return std::span<const std::byte>{std::data(source.Current.Buffer) + extent.Offset,
                                          extent.Length};
    }
};

MergeReadAhead::MergeReadAhead(std::span<const std::filesystem::path> inputs,
                               std::size_t decodeWorkers, ByteLimit budget, std::size_t batchBytes)
    : impl_{std::make_unique<Impl>(inputs, decodeWorkers, budget, batchBytes)}
{
}

MergeReadAhead::~MergeReadAhead() = default;

std::size_t MergeReadAhead::NumSources() const noexcept { return std::size(impl_->Sources); }

std::optional<std::span<const std::byte>> MergeReadAhead::Next(std::size_t sourceIndex)
{
    return impl_->Next(sourceIndex);
}

MergeReadAheadStats MergeReadAhead::Stats() const noexcept
{
    return MergeReadAheadStats{
        .ConsumerInputWaitNs = impl_->ConsumerInputWaitNs.load(std::memory_order_relaxed),
        .ProducerBudgetWaitNs = impl_->ProducerBudgetWaitNs.load(std::memory_order_relaxed),
        .IoReadNs = impl_->IoReadNs.load(std::memory_order_relaxed),
        .DecompressNs = impl_->DecompressNs.load(std::memory_order_relaxed),
        .PeakInFlightBytes = impl_->PeakInFlightBytes.load(std::memory_order_relaxed),
        .BudgetBytes = impl_->BudgetBytes,
    };
}

}  // namespace Samoa
}  // namespace PacBio
