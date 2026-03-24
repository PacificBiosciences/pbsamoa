#include <pbsamoa/io/BamZmwReader.hpp>

#include "ZmwUtils.hpp"

#include <pbsamoa/core/BamRecord.hpp>
#include <pbsamoa/core/SamHeader.hpp>
#include <pbsamoa/core/Tags.hpp>
#include <pbsamoa/index/ZmwIndex.hpp>
#include <pbsamoa/io/BamRecordReader.hpp>

#include <algorithm>
#include <condition_variable>
#include <deque>
#include <exception>
#include <expected>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {
namespace {

ZmwIdentity RecordZmwIdentity(const BamRecord& record)
{
    std::int32_t rgId{0};
    if (const TagValue* rgValue = record.Tags().Get(RG_TAG);
        (rgValue != nullptr) && std::holds_alternative<std::string>(*rgValue)) {
        rgId = ParseReadGroupId(std::get<std::string>(*rgValue));
    }
    return ZmwIdentity{rgId, ParseZmwFromName(record.Name())};
}

}  // namespace

struct BamZmwReader::Impl
{
    using QueueItem = std::expected<ZmwGroup, std::exception_ptr>;

    BamRecordReader reader_;
    SamHeader header_;
    ZmwIdentity currentZmw_{};
    std::optional<BamRecord> pending_;
    std::deque<QueueItem> queue_;
    std::mutex mutex_;
    std::condition_variable readyCv_;
    std::condition_variable spaceCv_;
    std::jthread producer_;
    std::size_t capacity_{0};
    bool done_{false};
    ZmwReaderMetrics metrics_{};

    explicit Impl(BamRecordReader reader, BamZmwReaderConfig config)
        : reader_{std::move(reader)}, header_{reader_.Header()}
    {
        if (config.PrefetchCapacityZmws == 0) {
            throw std::invalid_argument{"BamZmwReader prefetch capacity must be at least 1"};
        }
        capacity_ = config.PrefetchCapacityZmws;
        metrics_.ConfiguredCapacity = config.PrefetchCapacityZmws;
        producer_ = std::jthread{[this](std::stop_token stopToken) { ProducerLoop(stopToken); }};
    }

    ~Impl()
    {
        producer_.request_stop();
        readyCv_.notify_all();
        spaceCv_.notify_all();
        producer_.join();
    }

    Impl(const Impl&) = delete;
    Impl& operator=(const Impl&) = delete;
    Impl(Impl&&) = delete;
    Impl& operator=(Impl&&) = delete;

    bool EnsurePendingRecord()
    {
        if (!pending_) {
            pending_ = reader_.ReadRecord();
        }
        return static_cast<bool>(pending_);
    }

    std::optional<ZmwGroup> ReadNextGroup()
    {
        if (!EnsurePendingRecord()) {
            return std::nullopt;
        }

        ZmwGroup group{};
        group.zmw = RecordZmwIdentity(*pending_);
        group.records.push_back(std::move(*pending_));
        pending_.reset();

        while (auto rec{reader_.ReadRecord()}) {
            if (RecordZmwIdentity(*rec) == group.zmw) {
                group.records.push_back(std::move(*rec));
            } else {
                pending_ = std::move(rec);
                break;
            }
        }

        return group;
    }

    bool WaitForSpace(std::unique_lock<std::mutex>& lock, std::stop_token stopToken)
    {
        if (queue_.size() < capacity_) {
            return true;
        }

        metrics_.ProducerStalls++;
        spaceCv_.wait(lock, [this, stopToken] {
            return stopToken.stop_requested() || queue_.size() < capacity_;
        });
        return !stopToken.stop_requested();
    }

    bool PushItem(QueueItem item, std::stop_token stopToken)
    {
        std::unique_lock lock{mutex_};
        if (!WaitForSpace(lock, stopToken)) {
            done_ = true;
            readyCv_.notify_all();
            return false;
        }

        if (item.has_value()) {
            metrics_.GroupsProduced++;
        }
        queue_.push_back(std::move(item));
        metrics_.QueueDepth = queue_.size();
        metrics_.PeakQueueDepth = std::max(metrics_.PeakQueueDepth, metrics_.QueueDepth);
        readyCv_.notify_one();
        return true;
    }

    void FinishProduction()
    {
        {
            const std::lock_guard lock{mutex_};
            done_ = true;
        }
        readyCv_.notify_all();
    }

    void ProducerLoop(std::stop_token stopToken)
    {
        try {
            while (!stopToken.stop_requested()) {
                auto group = ReadNextGroup();
                if (!group) {
                    PushItem(QueueItem{std::unexpected{std::exception_ptr{}}}, stopToken);
                    FinishProduction();
                    return;
                }
                if (!PushItem(QueueItem{std::move(*group)}, stopToken)) {
                    return;
                }
            }
        } catch (...) {
            PushItem(QueueItem{std::unexpected{std::current_exception()}}, stopToken);
        }
        FinishProduction();
    }
};

BamZmwReader::BamZmwReader(BamRecordReader reader, BamZmwReaderConfig config)
    : impl_{std::make_unique<Impl>(std::move(reader), std::move(config))}
{
}

BamZmwReader::BamZmwReader(BamCollection collection, BamZmwReaderConfig config)
    : BamZmwReader{BamRecordReader{std::move(collection), config.Reader}, std::move(config)}
{
}

BamZmwReader::BamZmwReader(std::vector<std::filesystem::path> paths, BamZmwReaderConfig config)
    : BamZmwReader{BamRecordReader{std::move(paths), config.Reader}, std::move(config)}
{
}

BamZmwReader::BamZmwReader(std::vector<BamFile> files, BamZmwReaderConfig config)
    : BamZmwReader{BamRecordReader{std::move(files), config.Reader}, std::move(config)}
{
}

BamZmwReader::~BamZmwReader() = default;
BamZmwReader::BamZmwReader(BamZmwReader&&) noexcept = default;
BamZmwReader& BamZmwReader::operator=(BamZmwReader&&) noexcept = default;

bool BamZmwReader::GetNext(std::vector<BamRecord>& records)
{
    records.clear();

    std::unique_lock lock{impl_->mutex_};
    if (impl_->queue_.empty() && !impl_->done_) {
        impl_->metrics_.ConsumerStalls++;
        impl_->readyCv_.wait(lock, [this] { return !impl_->queue_.empty() || impl_->done_; });
    }

    if (impl_->queue_.empty()) {
        return false;
    }

    Impl::QueueItem item{std::move(impl_->queue_.front())};
    impl_->queue_.pop_front();
    impl_->metrics_.QueueDepth = impl_->queue_.size();
    if (item.has_value()) {
        impl_->metrics_.GroupsConsumed++;
        impl_->currentZmw_ = item->zmw;
    }
    lock.unlock();
    impl_->spaceCv_.notify_one();

    if (!item.has_value()) {
        if (const std::exception_ptr& ep = item.error(); ep) {
            std::rethrow_exception(ep);
        }
        return false;
    }

    records = std::move(item->records);
    return true;
}

ZmwIdentity BamZmwReader::CurrentZmw() const
{
    const std::lock_guard lock{impl_->mutex_};
    return impl_->currentZmw_;
}

const SamHeader& BamZmwReader::Header() const { return impl_->header_; }

std::int32_t BamZmwReader::NumZmws() const { return impl_->reader_.NumZmws(); }

ZmwReaderMetrics BamZmwReader::GetMetrics() const
{
    const std::lock_guard lock{impl_->mutex_};
    ZmwReaderMetrics metrics{impl_->metrics_};
    metrics.QueueDepth = impl_->queue_.size();
    metrics.Reader = impl_->reader_.GetMetrics();
    return metrics;
}

}  // namespace Samoa
}  // namespace PacBio
