#include <pbsamoa/io/BamZmwReader.hpp>

#include "ZmwUtils.hpp"

#include <pbsamoa/core/BamRecord.hpp>
#include <pbsamoa/core/SamHeader.hpp>
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
    std::size_t capacity_;
    bool done_{false};
    ZmwReaderMetrics metrics_{};

    static void RunProducerLoop(std::stop_token stopToken, Impl* self);

    explicit Impl(BamRecordReader reader, BamZmwReaderConfig config)
        : reader_{std::move(reader)}
        , header_{reader_.Header()}
        , capacity_{config.PrefetchCapacityZmws}
    {
        if (capacity_ == 0) {
            throw std::invalid_argument{"BamZmwReader prefetch capacity must be at least 1"};
        }
        metrics_.ConfiguredCapacity = capacity_;
        producer_ = std::jthread{&Impl::RunProducerLoop, this};
    }

    ~Impl()
    {
        producer_.request_stop();
        {
            const std::lock_guard lock{mutex_};
        }
        readyCv_.notify_all();
        spaceCv_.notify_all();
    }

    Impl(const Impl&) = delete;
    Impl& operator=(const Impl&) = delete;
    Impl(Impl&&) = delete;
    Impl& operator=(Impl&&) = delete;

    std::optional<ZmwGroup> ReadNextGroup()
    {
        if (!pending_) {
            pending_ = reader_.ReadRecord();
        }
        if (!pending_) {
            return std::nullopt;
        }

        ZmwGroup group{};
        group.zmw = ParseZmwIdentity(pending_->Name(), pending_->Tags());
        group.records.push_back(std::move(*pending_));
        pending_.reset();

        while (std::optional<BamRecord> rec{reader_.ReadRecord()}) {
            const ZmwIdentity recZmw{ParseZmwIdentity(rec->Name(), rec->Tags())};
            if (recZmw == group.zmw) {
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
        while ((queue_.size() >= capacity_) && !stopToken.stop_requested()) {
            spaceCv_.wait(lock);
        }
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

        if (item) {
            metrics_.GroupsProduced++;
        }
        queue_.push_back(std::move(item));
        metrics_.QueueDepth = queue_.size();
        metrics_.PeakQueueDepth = std::max(metrics_.PeakQueueDepth, metrics_.QueueDepth);
        readyCv_.notify_one();
        return true;
    }

    QueueItem PopFrontItem()
    {
        QueueItem item{std::move(queue_.front())};
        queue_.pop_front();
        metrics_.QueueDepth = queue_.size();
        return item;
    }

    void FinishProduction()
    {
        {
            const std::lock_guard lock{mutex_};
            done_ = true;
        }
        readyCv_.notify_all();
    }

    void PublishTerminalItem(std::exception_ptr error, std::stop_token stopToken)
    {
        PushItem(QueueItem{std::unexpected{error}}, stopToken);
        FinishProduction();
    }

    void ProducerLoop(std::stop_token stopToken)
    {
        try {
            while (!stopToken.stop_requested()) {
                std::optional<ZmwGroup> group{ReadNextGroup()};
                if (!group) {
                    PublishTerminalItem(std::exception_ptr{}, stopToken);
                    return;
                }
                if (!PushItem(QueueItem{std::move(*group)}, stopToken)) {
                    return;
                }
            }
        } catch (...) {
            PublishTerminalItem(std::current_exception(), stopToken);
            return;
        }
        FinishProduction();
    }
};

void BamZmwReader::Impl::RunProducerLoop(std::stop_token stopToken, Impl* self)
{
    self->ProducerLoop(stopToken);
}

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
    while (impl_->queue_.empty() && !impl_->done_) {
        impl_->metrics_.ConsumerStalls++;
        impl_->readyCv_.wait(lock);
    }

    if (impl_->queue_.empty()) {
        return false;
    }

    Impl::QueueItem item{impl_->PopFrontItem()};
    if (item) {
        impl_->metrics_.GroupsConsumed++;
        impl_->currentZmw_ = item->zmw;
    }
    lock.unlock();
    impl_->spaceCv_.notify_one();

    if (!item) {
        if (const std::exception_ptr error{item.error()}; error) {
            std::rethrow_exception(error);
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
