#include <pbsamoa/io/ZmiBamWriter.hpp>

#include <pbsamoa/io/BamWriter.hpp>
#include <pbsamoa/io/ZmiWriter.hpp>
#include "ZmwUtils.hpp"

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {

struct ZmiBamWriter::Impl
{
    ZmiWriter zmi;
    BamWriter bam;
    std::int32_t pendingRgId_{0};
    std::int32_t pendingZmw_{0};
    bool closed{false};

    Impl(const std::filesystem::path& bamPath, const SamHeader& header, int compressionLevel)
        : zmi{std::filesystem::path{bamPath.string() + ".zmi"}}
        , bam{bamPath, header, compressionLevel,
              [this](std::int64_t virtualOffset, std::span<const std::byte> /*rawData*/) {
                  zmi.AddRecord(pendingRgId_, pendingZmw_, virtualOffset);
              }}
    {
    }

    ~Impl() = default;
    Impl(const Impl&) = delete;
    Impl& operator=(const Impl&) = delete;
};

ZmiBamWriter::ZmiBamWriter(const std::filesystem::path& bamPath, const SamHeader& header,
                           int compressionLevel)
    : impl_{std::make_unique<Impl>(bamPath, header, compressionLevel)}
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

void ZmiBamWriter::Write(const BamRecord& record)
{
    impl_->pendingRgId_ = 0;
    impl_->pendingZmw_ = ParseZmwFromName(record.Name());
    impl_->bam.Write(record);
}

void ZmiBamWriter::Write(const RawRecord& view)
{
    impl_->pendingRgId_ = 0;
    impl_->pendingZmw_ = ParseZmwFromName(view.Name());
    impl_->bam.Write(view);
}

void ZmiBamWriter::WriteBatch(const RawRecordBatch& batch)
{
    for (std::size_t i{0}; i < batch.RecordCount(); ++i) {
        const RawRecord view{batch.RecordData(i)};
        Write(view);
    }
}

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
