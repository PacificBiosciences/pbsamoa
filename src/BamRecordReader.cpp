#include <pbsamoa/io/BamRecordReader.hpp>

#include <pbsamoa/core/BamRecord.hpp>
#include <pbsamoa/core/RawRecord.hpp>
#include <pbsamoa/core/Tags.hpp>
#include <pbsamoa/io/BamCollection.hpp>
#include <pbsamoa/io/BamRawReader.hpp>

#include <parallel/ThreadPool.h>
#include <rigtorp/SPSCQueue.hpp>

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <expected>
#include <filesystem>
#include <format>
#include <memory>
#include <mutex>
#include <optional>
#include <queue>
#include <thread>
#include <tuple>
#include <utility>
#include <variant>
#include <vector>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {
namespace {

using QueueItem = std::expected<BamRecord, std::exception_ptr>;
using DecodePool = PacBio::Parallel::ThreadPool<>;

std::shared_ptr<DecodePool> CreateDecodePool(std::size_t decodeWorkers)
{
    if (decodeWorkers == 0) {
        return {};
    }

    return std::make_shared<DecodePool>(DecodePool::Config{
        .NumThreads = decodeWorkers,
        .EnableMetrics = true,
    });
}

BamRawReader CreateViewReaderForCollection(const BamCollection& collection,
                                           BamRawReaderConfig config)
{
    // Keep source_ intact so Query() can build fresh readers later; the
    // sequential decode path consumes its own BamCollection copy.
    return BamRawReader{collection, std::move(config)};
}

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
    std::variant<std::filesystem::path, BamCollection> source_;
    BamRawReader viewReader_;
    SamHeader header_;
    std::shared_ptr<DecodePool> pool_;
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
        : source_{path}
        , viewReader_{path, config.RawReaderConfig}
        , header_{viewReader_.Header()}
        , pool_{CreateDecodePool(config.DecodeWorkers)}
        , queue_{config.OutputCapacity}
        , tagFilter_{config.TagFilter}
        , batchBudget_{config.BatchBudget}
        , parallelBgzf_{config.RawReaderConfig.BgzfWorkers > 0}
        , parallelDecode_{config.DecodeWorkers > 0}
    {
        StartProducer();
    }

    Impl(BamCollection collection, BamRecordReaderConfig config)
        : source_{std::move(collection)}
        , viewReader_{CreateViewReaderForCollection(std::get<BamCollection>(source_),
                                                    config.RawReaderConfig)}
        , header_{viewReader_.Header()}
        , pool_{CreateDecodePool(config.DecodeWorkers)}
        , queue_{config.OutputCapacity}
        , tagFilter_{config.TagFilter}
        , batchBudget_{config.BatchBudget}
        , parallelBgzf_{config.RawReaderConfig.BgzfWorkers > 0}
        , parallelDecode_{config.DecodeWorkers > 0}
    {
        StartProducer();
    }

    ~Impl() = default;
    Impl(const Impl&) = delete;
    Impl& operator=(const Impl&) = delete;
    Impl(Impl&&) = delete;
    Impl& operator=(Impl&&) = delete;

    void SignalConsumer()
    {
        // Pair the notify with a mutex acquisition so the consumer's wait/predicate
        // observes producer progress through the same synchronization point.
        {
            const std::lock_guard lock{readyMutex_};
        }
        readyCv_.notify_one();
    }

    void StartProducer()
    {
        producer_ = std::jthread{[this](std::stop_token stopToken) { ProducerLoop(stopToken); }};
    }

    QueueItem* WaitForQueueFront()
    {
        QueueItem* item{queue_.front()};
        if (item != nullptr) {
            return item;
        }

        counters_.consumerStalls.fetch_add(1, std::memory_order_relaxed);
        std::unique_lock lock{readyMutex_};
        readyCv_.wait(lock, [this]() {
            return queue_.front() != nullptr || done_.load(std::memory_order_acquire);
        });
        return queue_.front();
    }

    void ProducerLoop(std::stop_token stopToken)
    {
        try {
            while (!stopToken.stop_requested()) {
                const auto batchReadStart{std::chrono::steady_clock::now()};
                std::optional<RawRecordBatch> batch{viewReader_.ReadBatch(batchBudget_)};
                counters_.batchReadNs.fetch_add(
                    static_cast<std::uint64_t>(
                        std::chrono::duration_cast<std::chrono::nanoseconds>(
                            std::chrono::steady_clock::now() - batchReadStart)
                            .count()),
                    std::memory_order_relaxed);

                if (!batch) {
                    break;
                }

                const std::int32_t n{static_cast<std::int32_t>(batch->RecordCount())};
                std::vector<BamRecord> owned(n);

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

            done_.store(true, std::memory_order_release);
            while (!queue_.try_push(QueueItem{std::unexpected{std::exception_ptr{}}})) {
                if (stopToken.stop_requested()) {
                    return;
                }
                std::this_thread::yield();
            }
            SignalConsumer();
        } catch (...) {
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

struct BamRecordReader::QueryRange::Impl
{
    struct Source
    {
        std::unique_ptr<BamRawReader> reader;
        std::optional<BamRawReader::QueryRange> range;
        std::optional<BamRawReader::QueryRange::Iterator> iter;
        std::optional<BamRawReader::QueryRange::Iterator> end;
        std::optional<RawRecord> current;
        std::size_t order{0};
    };

    struct SourceCompare
    {
        bool operator()(const Source* lhs, const Source* rhs) const
        {
            const auto lhsKey = std::tuple{lhs->current->RefId(), lhs->current->Pos(),
                                           lhs->current->Name(), lhs->order};
            const auto rhsKey = std::tuple{rhs->current->RefId(), rhs->current->Pos(),
                                           rhs->current->Name(), rhs->order};
            return lhsKey > rhsKey;
        }
    };

    std::vector<Source> sources_;
    std::priority_queue<Source*, std::vector<Source*>, SourceCompare> ready_;
    std::variant<std::monostate, DropTags, KeepTags> tagFilter_;

    Impl(BamCollection collection, const GenomicInterval& interval,
         std::variant<std::monostate, DropTags, KeepTags> tagFilter)
        : tagFilter_{std::move(tagFilter)}
    {
        const std::int32_t refId{collection.Header().ReferenceId(interval.Name())};
        if (refId < 0) {
            throw std::runtime_error{
                std::format("reference '{}' not found in BAM header", interval.Name())};
        }

        sources_.reserve(collection.Size());
        for (std::size_t i{0}; i < collection.Size(); ++i) {
            const BamFile& file{collection.Files()[i]};
            const std::filesystem::path baiPath{file.StandardIndexFilename()};
            if (!std::filesystem::exists(baiPath)) {
                throw std::runtime_error{std::format("index file not found: {}", baiPath.string())};
            }

            const BaiIndex index{BaiIndex::FromFile(baiPath)};
            sources_.push_back(Source{
                .reader = std::make_unique<BamRawReader>(file.Filename()),
                .range = std::nullopt,
                .iter = std::nullopt,
                .end = std::nullopt,
                .current = std::nullopt,
                .order = i,
            });
            Source& source{sources_.back()};
            source.range.emplace(source.reader.get(),
                                 index.Query(refId, interval.Start(), interval.Stop()), refId,
                                 interval.Start(), interval.Stop());
            source.iter.emplace(source.range->begin());
            source.end.emplace(source.range->end());
            if ((*source.iter) == (*source.end)) {
                continue;
            }
            source.current.emplace(**source.iter);
            ready_.push(&source);
        }
    }
};

BamRecordReader::BamRecordReader(const std::filesystem::path& path, BamRecordReaderConfig config)
    : impl_{std::make_unique<Impl>(path, std::move(config))}
{
}

BamRecordReader::BamRecordReader(BamCollection collection, BamRecordReaderConfig config)
    : impl_{std::make_unique<Impl>(std::move(collection), std::move(config))}
{
}

BamRecordReader::BamRecordReader(std::vector<std::filesystem::path> paths,
                                 BamRecordReaderConfig config)
    : BamRecordReader{BamCollection{std::move(paths)}, std::move(config)}
{
}

BamRecordReader::BamRecordReader(std::vector<BamFile> files, BamRecordReaderConfig config)
    : BamRecordReader{BamCollection{std::move(files)}, std::move(config)}
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

    QueueItem* item{impl_->WaitForQueueFront()};
    if (item == nullptr) {
        impl_->eof_ = true;
        return std::nullopt;
    }

    if (!item->has_value()) {
        const std::exception_ptr& ep{item->error()};
        impl_->queue_.pop();
        if (ep) {
            std::rethrow_exception(ep);
        }
        impl_->eof_ = true;
        return std::nullopt;
    }

    BamRecord record{std::move(item->value())};
    impl_->queue_.pop();
    impl_->counters_.recordsConsumed.fetch_add(1, std::memory_order_relaxed);
    return record;
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
    if (!current_) {
        reader_ = nullptr;
    }
}

const BamRecord& BamRecordReader::RecordRange::Iterator::operator*() const { return *current_; }

const BamRecord* BamRecordReader::RecordRange::Iterator::operator->() const { return &*current_; }

BamRecordReader::RecordRange::Iterator& BamRecordReader::RecordRange::Iterator::operator++()
{
    current_ = reader_->ReadRecord();
    if (!current_) {
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

// --- QueryRange ---

BamRecordReader::QueryRange::QueryRange(std::filesystem::path path, GenomicInterval interval,
                                        std::variant<std::monostate, DropTags, KeepTags> tagFilter)
    : QueryRange{BamCollection{std::move(path)}, std::move(interval), std::move(tagFilter)}
{
}

BamRecordReader::QueryRange::QueryRange(BamCollection collection, GenomicInterval interval,
                                        std::variant<std::monostate, DropTags, KeepTags> tagFilter)
    : impl_{std::make_unique<Impl>(std::move(collection), interval, std::move(tagFilter))}
{
}

BamRecordReader::QueryRange::~QueryRange() = default;
BamRecordReader::QueryRange::QueryRange(QueryRange&&) noexcept = default;
BamRecordReader::QueryRange& BamRecordReader::QueryRange::operator=(QueryRange&&) noexcept =
    default;

BamRecordReader::QueryRange::Iterator::Iterator() = default;

BamRecordReader::QueryRange::Iterator::Iterator(QueryRange* range) : range_{range}
{
    current_ = range_->ReadRecord();
    if (!current_) {
        range_ = nullptr;
    }
}

const BamRecord& BamRecordReader::QueryRange::Iterator::operator*() const { return *current_; }

const BamRecord* BamRecordReader::QueryRange::Iterator::operator->() const { return &*current_; }

BamRecordReader::QueryRange::Iterator& BamRecordReader::QueryRange::Iterator::operator++()
{
    current_ = range_->ReadRecord();
    if (!current_) {
        range_ = nullptr;
    }
    return *this;
}

void BamRecordReader::QueryRange::Iterator::operator++(int) { ++(*this); }

bool BamRecordReader::QueryRange::Iterator::operator==(const Iterator& other) const
{
    return range_ == other.range_;
}

BamRecordReader::QueryRange::Iterator BamRecordReader::QueryRange::begin()
{
    return Iterator{this};
}

BamRecordReader::QueryRange::Iterator BamRecordReader::QueryRange::end() { return Iterator{}; }

std::optional<BamRecord> BamRecordReader::QueryRange::ReadRecord()
{
    if (impl_->ready_.empty()) {
        return std::nullopt;
    }

    QueryRange::Impl::Source* source{impl_->ready_.top()};
    impl_->ready_.pop();

    BamRecord record{DecodeView(*source->current, impl_->tagFilter_)};
    ++(*source->iter);
    if ((*source->iter) != (*source->end)) {
        source->current.emplace(**source->iter);
        impl_->ready_.push(source);
    } else {
        source->current.reset();
    }
    return record;
}

BamRecordReader::QueryRange BamRecordReader::Query(std::string_view refName, std::int32_t beg,
                                                   std::int32_t end)
{
    const GenomicInterval interval{std::string{refName}, beg, end};
    return Query(interval);
}

BamRecordReader::QueryRange BamRecordReader::Query(const GenomicInterval& interval)
{
    return std::visit(
        [this, &interval](const auto& source) -> QueryRange {
            return QueryRange{source, interval, impl_->tagFilter_};
        },
        impl_->source_);
}

ReaderMetrics BamRecordReader::GetMetrics() const
{
    ReaderMetrics m{};

    m.Bgzf = impl_->viewReader_.GetMetrics();

    if (impl_->pool_) {
        const auto poolSnap{impl_->pool_->GetMetrics()};
        m.Decode.Pool.QueueDepth = poolSnap.CurrentQueueDepth;
        m.Decode.Pool.PeakQueueDepth = poolSnap.PeakQueueDepth;
        m.Decode.Pool.ActiveWorkers = poolSnap.CurrentActiveTasks;
        m.Decode.Pool.PeakActiveWorkers = poolSnap.PeakActiveTasks;
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
