#include <pbsamoa/io/BamZmwReader.hpp>

#include "ZmwUtils.hpp"

#include <pbsamoa/core/BamRecord.hpp>
#include <pbsamoa/core/SamHeader.hpp>
#include <pbsamoa/core/Tags.hpp>
#include <pbsamoa/index/ZmwIndex.hpp>
#include <pbsamoa/io/BamRecordReader.hpp>

#include <optional>
#include <string>
#include <utility>
#include <vector>

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

    bool EnsurePendingRecord()
    {
        if (!pending_) {
            pending_ = reader_.ReadRecord();
        }
        return static_cast<bool>(pending_);
    }
};

BamZmwReader::BamZmwReader(BamRecordReader reader)
    : impl_{std::make_unique<Impl>(std::move(reader))}
{
}

BamZmwReader::BamZmwReader(BamCollection collection, BamRecordReaderConfig config)
    : BamZmwReader{BamRecordReader{std::move(collection), std::move(config)}}
{
}

BamZmwReader::BamZmwReader(std::vector<std::filesystem::path> paths, BamRecordReaderConfig config)
    : BamZmwReader{BamRecordReader{std::move(paths), std::move(config)}}
{
}

BamZmwReader::BamZmwReader(std::vector<BamFile> files, BamRecordReaderConfig config)
    : BamZmwReader{BamRecordReader{std::move(files), std::move(config)}}
{
}

BamZmwReader::~BamZmwReader() = default;
BamZmwReader::BamZmwReader(BamZmwReader&&) noexcept = default;
BamZmwReader& BamZmwReader::operator=(BamZmwReader&&) noexcept = default;

bool BamZmwReader::GetNext(std::vector<BamRecord>& records)
{
    records.clear();

    if (!impl_->EnsurePendingRecord()) {
        return false;
    }

    const ZmwIdentity groupIdentity{RecordZmwIdentity(*impl_->pending_)};
    impl_->currentZmw_ = groupIdentity;
    records.push_back(std::move(*impl_->pending_));
    impl_->pending_.reset();

    while (auto rec{impl_->reader_.ReadRecord()}) {
        if (RecordZmwIdentity(*rec) == groupIdentity) {
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
