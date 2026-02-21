#include <pbsamoa/io/ZmiBamWriter.hpp>

#include "CramInternal.hpp"
#include "ZmwUtils.hpp"

#include <pbsamoa/core/RawRecord.hpp>
#include <pbsamoa/io/BamWriter.hpp>
#include <pbsamoa/io/ZmiWriter.hpp>

#include <string>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {
namespace {

std::int32_t ReadGroupIdFromTags(const TagMap& tags)
{
    const TagValue* rgValue{tags.Get(RG_TAG)};
    if ((rgValue == nullptr) || !std::holds_alternative<std::string>(*rgValue)) {
        return 0;
    }
    return ParseReadGroupId(std::get<std::string>(*rgValue));
}

void AddRawRecordToIndex(ZmiWriter& zmi, std::int64_t virtualOffset,
                         std::span<const std::byte> rawData)
{
    const RawRecord rawRecord{rawData};
    const TagMap tags{rawRecord.ParseTags()};
    zmi.AddRecord(ReadGroupIdFromTags(tags), ParseZmwFromName(rawRecord.Name()), virtualOffset);
}

}  // namespace

struct ZmiBamWriter::Impl
{
    ZmiWriter zmi;
    BamWriter bam;
    bool closed{false};

    Impl(const std::filesystem::path& bamPath, const SamHeader& header, const BamWriterConfig& cfg)
        : zmi{std::filesystem::path{bamPath.string() + ".zmi"},
              ZmiWriterConfig{.UseTempFile = cfg.UseTempFile}}
        , bam{bamPath, header, cfg,
              [this](std::int64_t virtualOffset, std::span<const std::byte> rawData) {
                  AddRawRecordToIndex(zmi, virtualOffset, rawData);
              }}
    {
    }

    ~Impl() = default;
    Impl(const Impl&) = delete;
    Impl& operator=(const Impl&) = delete;
};

ZmiBamWriter::ZmiBamWriter(const std::filesystem::path& bamPath, const SamHeader& header,
                           const BamWriterConfig& config)
    : impl_{std::make_unique<Impl>(bamPath, header, config)}
{
}

ZmiBamWriter::~ZmiBamWriter()
{
    if ((impl_ != nullptr) && (!impl_->closed)) {
        Close();
    }
}

ZmiBamWriter::ZmiBamWriter(ZmiBamWriter&&) noexcept = default;
ZmiBamWriter& ZmiBamWriter::operator=(ZmiBamWriter&&) noexcept = default;

void ZmiBamWriter::Write(const BamRecord& record) { impl_->bam.Write(record); }

void ZmiBamWriter::Write(const RawRecord& view) { impl_->bam.Write(view); }

void ZmiBamWriter::WriteBatch(const RawRecordBatch& batch) { impl_->bam.WriteBatch(batch); }

void ZmiBamWriter::Close()
{
    if (impl_->closed) {
        return;
    }
    // BamWriter must close first (flushes remaining records + fires callbacks)
    impl_->bam.Close();
    // Then ZmiWriter closes (finalizes the index)
    impl_->zmi.Close();
    impl_->closed = true;
}

}  // namespace Samoa
}  // namespace PacBio
