#include <pbsamoa/io/ZmiBamWriter.hpp>

#include "PathUtils.hpp"
#include "ZmwUtils.hpp"

#include <pbsamoa/core/RawRecord.hpp>
#include <pbsamoa/io/BamWriter.hpp>
#include <pbsamoa/io/ZmiWriter.hpp>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {

struct ZmiBamWriter::Impl
{
    ZmiWriter zmi;
    BamWriter bam;
    bool closed{false};

    Impl(const std::filesystem::path& bamPath, const SamHeader& header, const BamWriterConfig& cfg)
        : zmi{SidecarPath(bamPath, ".zmi"), ZmiWriterConfig{.UseTempFile = cfg.UseTempFile}}
        , bam{bamPath, header, cfg,
              [this](std::int64_t virtualOffset, std::span<const std::byte> rawData) {
                  const RawRecord rawRecord{rawData};
                  const ZmwIdentity identity{
                      ParseZmwIdentity(rawRecord.Name(), rawRecord.ParseTags())};
                  zmi.AddRecord(identity.rgId, identity.zmw, virtualOffset);
              }}
    {
    }

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
    if (impl_ && !impl_->closed) {
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
