#include <pbsamoa/io/BamZmwReader.hpp>

#include "ZmwUtils.hpp"

#include <pbsamoa/core/BamRecord.hpp>
#include <pbsamoa/core/SamHeader.hpp>
#include <pbsamoa/index/ZmwIndex.hpp>
#include <pbsamoa/io/BamRecordReader.hpp>

#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <cstdint>

namespace PacBio {
namespace Samoa {

struct BamZmwReader::Impl
{
    BamRecordReader reader_;
    SamHeader header_;
    ZmwIdentity currentZmw_{};
    std::optional<BamRecord> pending_;

    explicit Impl(BamRecordReader reader) : reader_{std::move(reader)}, header_{reader_.Header()} {}

    ~Impl() = default;
    Impl(const Impl&) = delete;
    Impl& operator=(const Impl&) = delete;
    Impl(Impl&&) = delete;
    Impl& operator=(Impl&&) = delete;
};

BamZmwReader::BamZmwReader(BamRecordReader reader)
    : impl_{std::make_unique<Impl>(std::move(reader))}
{
}

BamZmwReader::~BamZmwReader() = default;
BamZmwReader::BamZmwReader(BamZmwReader&&) noexcept = default;
BamZmwReader& BamZmwReader::operator=(BamZmwReader&&) noexcept = default;

bool BamZmwReader::GetNext(std::vector<BamRecord>& records)
{
    records.clear();

    // Use pending record from last call, or read new one
    if (!impl_->pending_.has_value()) {
        impl_->pending_ = impl_->reader_.ReadRecord();
    }
    if (!impl_->pending_.has_value()) {
        return false;
    }

    // Start group with pending record — use full movie/zmw prefix for grouping
    // to distinguish records from different movies with the same hole number
    const std::string groupKey{MovieZmwPrefix(impl_->pending_->Name())};
    const std::int32_t zmw{ParseZmwFromName(impl_->pending_->Name())};
    impl_->currentZmw_ = ZmwIdentity{0, zmw};
    records.push_back(std::move(*impl_->pending_));
    impl_->pending_.reset();

    // Accumulate consecutive records with same movie/zmw
    while (true) {
        auto rec{impl_->reader_.ReadRecord()};
        if (!rec.has_value()) {
            break;
        }
        if (MovieZmwPrefix(rec->Name()) == groupKey) {
            records.push_back(std::move(*rec));
        } else {
            impl_->pending_ = std::move(rec);
            break;
        }
    }

    return true;
}

ZmwIdentity BamZmwReader::CurrentZmw() const { return impl_->currentZmw_; }

const SamHeader& BamZmwReader::Header() const { return impl_->header_; }

}  // namespace Samoa
}  // namespace PacBio
