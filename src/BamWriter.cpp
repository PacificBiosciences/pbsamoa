#include <pbsamoa/io/BamWriter.hpp>

#include <pbsamoa/core/Bgzf.hpp>

#include <algorithm>
#include <array>
#include <format>
#include <functional>
#include <limits>
#include <stdexcept>
#include <vector>

#include <cstdint>

namespace PacBio {
namespace Samoa {

struct BamWriter::Impl
{
    BgzfWriter bgzf;
    IndexCallback callback{};
    bool closed{false};

    Impl(const std::filesystem::path& path, const SamHeader& header, int compressionLevel)
        : bgzf{path, compressionLevel}
    {
        const std::vector<std::byte> headerBlock{header.ToBamHeaderBlock()};
        bgzf.Write(headerBlock);
    }

    Impl(const std::filesystem::path& path, const SamHeader& header, int compressionLevel,
         IndexCallback cb)
        : bgzf{path, compressionLevel}, callback{std::move(cb)}
    {
        const std::vector<std::byte> headerBlock{header.ToBamHeaderBlock()};
        bgzf.Write(headerBlock);
    }
};

BamWriter::BamWriter(const std::filesystem::path& path, const SamHeader& header,
                     int compressionLevel)
    : impl_{std::make_unique<Impl>(path, header, compressionLevel)}
{
}

BamWriter::BamWriter(const std::filesystem::path& path, const SamHeader& header,
                     int compressionLevel, IndexCallback callback)
    : impl_{std::make_unique<Impl>(path, header, compressionLevel, std::move(callback))}
{
}

BamWriter::~BamWriter()
{
    if ((impl_ != nullptr) && (!impl_->closed)) {
        Close();
    }
}

BamWriter::BamWriter(BamWriter&&) noexcept = default;
BamWriter& BamWriter::operator=(BamWriter&&) noexcept = default;

void BamWriter::Write(const BamRecord& record)
{
    const VirtualOffset offset{impl_->bgzf.Tell()};
    const std::vector<std::byte> serialized{record.SerializeToBam()};
    if (std::size(serialized) >
        static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max())) {
        throw std::invalid_argument{std::format(
            "BamWriter::Write: serialized record size ({}) exceeds BAM block_size limit",
            std::size(serialized))};
    }
    const std::uint32_t blockSize{static_cast<std::uint32_t>(std::size(serialized))};

    std::array<std::byte, 4> prefix;
    std::ranges::copy_n(reinterpret_cast<const std::byte*>(&blockSize), sizeof(blockSize),
                        std::data(prefix));
    impl_->bgzf.Write(prefix);
    impl_->bgzf.Write(serialized);

    if (impl_->callback) {
        impl_->callback(offset.Value(), serialized);
    }
}

void BamWriter::Write(std::span<const std::byte> rawData)
{
    const VirtualOffset offset{impl_->bgzf.Tell()};
    if (std::size(rawData) > static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max())) {
        throw std::invalid_argument{
            std::format("BamWriter::Write: raw record size ({}) exceeds BAM block_size limit",
                        std::size(rawData))};
    }
    const std::uint32_t blockSize{static_cast<std::uint32_t>(std::size(rawData))};

    std::array<std::byte, 4> prefix;
    std::ranges::copy_n(reinterpret_cast<const std::byte*>(&blockSize), sizeof(blockSize),
                        std::data(prefix));
    impl_->bgzf.Write(prefix);
    impl_->bgzf.Write(rawData);

    if (impl_->callback) {
        impl_->callback(offset.Value(), rawData);
    }
}

void BamWriter::Write(const RawRecord& byteView) { Write(byteView.RawData()); }

void BamWriter::WriteBatch(const RawRecordBatch& batch)
{
    for (std::size_t i{0}; i < batch.RecordCount(); ++i) {
        Write(batch.RecordData(i));
    }
}

void BamWriter::Close()
{
    if (impl_->closed) {
        return;
    }
    impl_->bgzf.Close();
    impl_->closed = true;
}

}  // namespace Samoa
}  // namespace PacBio
