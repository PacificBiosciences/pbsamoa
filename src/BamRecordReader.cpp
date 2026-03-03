#include <pbsamoa/io/BamRecordReader.hpp>

#include <pbsamoa/core/BamRecord.hpp>
#include <pbsamoa/core/RawRecord.hpp>
#include <pbsamoa/core/Tags.hpp>
#include <pbsamoa/io/BamRawReader.hpp>

#include <parallel/ThreadPool.h>
#include <rigtorp/SPSCQueue.hpp>

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <expected>
#include <memory>
#include <mutex>
#include <optional>
#include <thread>
#include <utility>
#include <variant>
#include <vector>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {
namespace {

using QueueItem = std::expected<BamRecord, std::exception_ptr>;

BamRecord DecodeView(const RawRecord& view,
                     const std::variant<std::monostate, DropTags, KeepTags>& tagFilter)
{
    return std::visit(
        [&](const auto& filter) -> BamRecord {
            using T = std::remove_cvref_t<decltype(filter)>;
            if constexpr (std::is_same_v<T, std::monostate>) {
                return view.ToOwned();
            } else {
                return view.ToOwned(filter);
            }
        },
        tagFilter);
}

}  // namespace

/// \brief Always-on atomic counters for decode pipeline introspection.
struct DecodeCounters
{
    std::atomic<std::uint64_t> batchesDecoded{0};
    std::atomic<std::uint64_t> recordsDecoded{0};
    std::atomic<std::uint64_t> recordsProduced{0};
    std::atomic<std::uint64_t> recordsConsumed{0};
    std::atomic<std::uint64_t> producerStalls{0};
    std::atomic<std::uint64_t> consumerStalls{0};
    std::atomic<std::uint64_t> decodeNs{0};
    std::atomic<std::uint64_t> batchReadNs{0};
};

struct BamRecordReader::Impl
{
    BamRawReader viewReader_;
    SamHeader header_;
    std::shared_ptr<PacBio::Parallel::ThreadPool<>> pool_;
    rigtorp::SPSCQueue<QueueItem> queue_;
    std::variant<std::monostate, DropTags, KeepTags> tagFilter_;
    ByteLimit batchBudget_;
    bool eof_{false};
    bool parallelBgzf_;
    bool parallelDecode_;
    DecodeCounters counters_;
    std::mutex readyMutex_;
    std::condition_variable readyCv_;
    std::atomic<bool> done_{false};
    std::jthread producer_;

    Impl(const std::filesystem::path& path, BamRecordReaderConfig config)
        : viewReader_{path, config.RawReaderConfig}
        , header_{viewReader_.Header()}
        , pool_{config.DecodeWorkers > 0 ? std::make_shared<PacBio::Parallel::ThreadPool<>>(
                                               PacBio::Parallel::ThreadPool<>::Config{
                                                   .NumThreads = config.DecodeWorkers,
                                                   .EnableMetrics = true,
                                               })
                                         : std::shared_ptr<PacBio::Parallel::ThreadPool<>>{}}
        , queue_{config.OutputCapacity}
        , tagFilter_{config.TagFilter}
        , batchBudget_{config.BatchBudget}
        , parallelBgzf_{config.RawReaderConfig.BgzfWorkers > 0}
        , parallelDecode_{config.DecodeWorkers > 0}
    {
        producer_ = std::jthread{[this](std::stop_token stopToken) { ProducerLoop(stopToken); }};
    }

    ~Impl() = default;
    Impl(const Impl&) = delete;
    Impl& operator=(const Impl&) = delete;
    Impl(Impl&&) = delete;
    Impl& operator=(Impl&&) = delete;

    /// \brief Acquire-release the mutex so the consumer sees all prior writes,
    ///        then wake one waiter. The lock scope is intentionally empty: the
    ///        unlock acts as a release fence pairing with the consumer's acquire.
    void SignalConsumer()
    {
        {
            const std::lock_guard lock{readyMutex_};
        }
        readyCv_.notify_one();
    }

    void ProducerLoop(std::stop_token stopToken)
    {
        try {
            while (!stopToken.stop_requested()) {
                // Time the batch read from underlying reader
                const auto batchReadStart{std::chrono::steady_clock::now()};
                std::optional<RawRecordBatch> batch{viewReader_.ReadBatch(batchBudget_)};
                counters_.batchReadNs.fetch_add(
                    static_cast<std::uint64_t>(
                        std::chrono::duration_cast<std::chrono::nanoseconds>(
                            std::chrono::steady_clock::now() - batchReadStart)
                            .count()),
                    std::memory_order_relaxed);

                if (!batch.has_value()) {
                    break;
                }

                const std::int32_t n{static_cast<std::int32_t>(batch->RecordCount())};
                std::vector<BamRecord> owned(n);

                // Time the parallel decode
                const auto decodeStart{std::chrono::steady_clock::now()};
                PacBio::Parallel::Dispatch(
                    pool_,
                    [&](std::int32_t i) {
                        const RawRecord view{batch->RecordData(i)};
                        owned[i] = DecodeView(view, tagFilter_);
                    },
                    n);
                counters_.decodeNs.fetch_add(
                    static_cast<std::uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(
                                                   std::chrono::steady_clock::now() - decodeStart)
                                                   .count()),
                    std::memory_order_relaxed);

                counters_.batchesDecoded.fetch_add(1, std::memory_order_relaxed);
                counters_.recordsDecoded.fetch_add(static_cast<std::uint64_t>(n),
                                                   std::memory_order_relaxed);

                for (auto& rec : owned) {
                    if (stopToken.stop_requested()) {
                        return;
                    }
                    while (!queue_.try_push(QueueItem{std::move(rec)})) {
                        if (stopToken.stop_requested()) {
                            return;
                        }
                        counters_.producerStalls.fetch_add(1, std::memory_order_relaxed);
                        std::this_thread::yield();
                    }
                    counters_.recordsProduced.fetch_add(1, std::memory_order_relaxed);
                    SignalConsumer();
                }
            }

            // EOF sentinel
            done_.store(true, std::memory_order_release);
            while (!queue_.try_push(QueueItem{std::unexpected{std::exception_ptr{}}})) {
                if (stopToken.stop_requested()) {
                    return;
                }
                std::this_thread::yield();
            }
            SignalConsumer();
        } catch (...) {
            // Error sentinel
            done_.store(true, std::memory_order_release);
            while (!queue_.try_push(QueueItem{std::unexpected{std::current_exception()}})) {
                if (stopToken.stop_requested()) {
                    return;
                }
                std::this_thread::yield();
            }
            SignalConsumer();
        }
    }
};

BamRecordReader::BamRecordReader(const std::filesystem::path& path, BamRecordReaderConfig config)
    : impl_{std::make_unique<Impl>(path, std::move(config))}
{
}

BamRecordReader::~BamRecordReader() = default;
BamRecordReader::BamRecordReader(BamRecordReader&&) noexcept = default;
BamRecordReader& BamRecordReader::operator=(BamRecordReader&&) noexcept = default;

const SamHeader& BamRecordReader::Header() const { return impl_->header_; }

std::optional<BamRecord> BamRecordReader::ReadRecord()
{
    if (impl_->eof_) {
        return std::nullopt;
    }

    // Fast path: spin on SPSC queue
    QueueItem* result{impl_->queue_.front()};
    if (result == nullptr) {
        // Slow path: wait on condition variable
        impl_->counters_.consumerStalls.fetch_add(1, std::memory_order_relaxed);
        std::unique_lock lock{impl_->readyMutex_};
        impl_->readyCv_.wait(lock, [this]() {
            return impl_->queue_.front() != nullptr || impl_->done_.load(std::memory_order_acquire);
        });
        result = impl_->queue_.front();
        if (result == nullptr) {
            impl_->eof_ = true;
            return std::nullopt;
        }
    }

    if (result->has_value()) {
        BamRecord rec{std::move(result->value())};
        impl_->queue_.pop();
        impl_->counters_.recordsConsumed.fetch_add(1, std::memory_order_relaxed);
        return rec;
    }

    // Unexpected: EOF or error
    const std::exception_ptr& ep{result->error()};
    impl_->queue_.pop();
    if (ep) {
        std::rethrow_exception(ep);
    }
    impl_->eof_ = true;
    return std::nullopt;
}

// --- RecordRange ---

BamRecordReader::RecordRange::RecordRange(BamRecordReader* reader) : reader_{reader} {}

BamRecordReader::RecordRange::Iterator BamRecordReader::RecordRange::begin()
{
    return Iterator{reader_};
}

BamRecordReader::RecordRange::Iterator BamRecordReader::RecordRange::end() { return Iterator{}; }

// --- Iterator ---

BamRecordReader::RecordRange::Iterator::Iterator() = default;

BamRecordReader::RecordRange::Iterator::Iterator(BamRecordReader* reader) : reader_{reader}
{
    current_ = reader_->ReadRecord();
    if (!current_.has_value()) {
        reader_ = nullptr;
    }
}

const BamRecord& BamRecordReader::RecordRange::Iterator::operator*() const { return *current_; }

const BamRecord* BamRecordReader::RecordRange::Iterator::operator->() const { return &*current_; }

BamRecordReader::RecordRange::Iterator& BamRecordReader::RecordRange::Iterator::operator++()
{
    current_ = reader_->ReadRecord();
    if (!current_.has_value()) {
        reader_ = nullptr;
    }
    return *this;
}

void BamRecordReader::RecordRange::Iterator::operator++(int) { ++(*this); }

bool BamRecordReader::RecordRange::Iterator::operator==(const Iterator& other) const
{
    return reader_ == other.reader_;
}

BamRecordReader::RecordRange BamRecordReader::Records() { return RecordRange{this}; }

ReaderMetrics BamRecordReader::GetMetrics() const
{
    ReaderMetrics m{};

    // BGZF layer metrics from underlying raw reader
    m.Bgzf = impl_->viewReader_.GetMetrics();

    // Decode layer metrics
    if (impl_->pool_) {
        const auto poolSnap{impl_->pool_->GetMetrics()};
        m.Decode.PoolQueueDepth = poolSnap.CurrentQueueDepth;
        m.Decode.PoolPeakQueueDepth = poolSnap.PeakQueueDepth;
        m.Decode.PoolActiveWorkers = poolSnap.CurrentActiveTasks;
        m.Decode.PoolPeakActiveWorkers = poolSnap.PeakActiveTasks;
    }

    m.Decode.BatchesDecoded = impl_->counters_.batchesDecoded.load(std::memory_order_relaxed);
    m.Decode.RecordsDecoded = impl_->counters_.recordsDecoded.load(std::memory_order_relaxed);
    m.Decode.RecordsProduced = impl_->counters_.recordsProduced.load(std::memory_order_relaxed);
    m.Decode.RecordsConsumed = impl_->counters_.recordsConsumed.load(std::memory_order_relaxed);
    m.Decode.ProducerStalls = impl_->counters_.producerStalls.load(std::memory_order_relaxed);
    m.Decode.ConsumerStalls = impl_->counters_.consumerStalls.load(std::memory_order_relaxed);
    m.Decode.DecodeNs = impl_->counters_.decodeNs.load(std::memory_order_relaxed);
    m.Decode.BatchReadNs = impl_->counters_.batchReadNs.load(std::memory_order_relaxed);

    m.TotalRecordsRead = impl_->counters_.recordsConsumed.load(std::memory_order_relaxed);
    m.ParallelBgzf = impl_->parallelBgzf_;
    m.ParallelDecode = impl_->parallelDecode_;

    return m;
}

}  // namespace Samoa
}  // namespace PacBio
