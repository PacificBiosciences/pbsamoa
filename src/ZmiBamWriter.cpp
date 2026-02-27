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
    bool closed{false};

    static std::int32_t ParseZmwFromRawRecord(std::span<const std::byte> rawData)
    {
        // rawData starts at BAM refID; l_read_name is byte 8, read name starts at
        // byte 32.
        constexpr std::size_t FIXED_FIELDS_SIZE{32};
        constexpr std::size_t NAME_LENGTH_OFFSET{8};

        if (std::size(rawData) <= NAME_LENGTH_OFFSET) {
            return 0;
        }
        const std::uint8_t lReadName{std::to_integer<std::uint8_t>(rawData[NAME_LENGTH_OFFSET])};
        if ((lReadName == 0U) || (std::size(rawData) < (FIXED_FIELDS_SIZE + lReadName))) {
            return 0;
        }

        const char* namePtr{reinterpret_cast<const char*>(std::data(rawData) + FIXED_FIELDS_SIZE)};
        const std::string_view readName{namePtr, static_cast<std::size_t>(lReadName - 1U)};
        return ParseZmwFromName(readName);
    }

    Impl(const std::filesystem::path& bamPath, const SamHeader& header, const BamWriterConfig& cfg)
        : zmi{std::filesystem::path{bamPath.string() + ".zmi"}}
        , bam{bamPath, header, cfg,
              [this](std::int64_t virtualOffset, std::span<const std::byte> rawData) {
                  const std::int32_t zmw{ParseZmwFromRawRecord(rawData)};
                  zmi.AddRecord(0, zmw, virtualOffset);
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
