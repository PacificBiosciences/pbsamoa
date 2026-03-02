#include <pbsamoa/io/BamWriter.hpp>

#include <algorithm>
#include <array>
#include <atomic>
#include <format>
#include <functional>
#include <limits>
#include <stdexcept>
#include <vector>

#include <cstdint>

namespace PacBio {
namespace Samoa {

namespace {

std::vector<std::byte> BuildCombinedRecordPayload(std::span<const std::byte> rawData)
{
    if (std::size(rawData) > static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max())) {
        throw std::invalid_argument{
            std::format("BamWriter::Write: raw record size ({}) exceeds BAM block_size limit",
                        std::size(rawData))};
    }

    const std::uint32_t blockSize = std::size(rawData);

    std::vector<std::byte> combined{};
    combined.resize(4U + std::size(rawData));
    combined[0] = static_cast<std::byte>(blockSize & 0xFFU);
    combined[1] = static_cast<std::byte>((blockSize >> 8U) & 0xFFU);
    combined[2] = static_cast<std::byte>((blockSize >> 16U) & 0xFFU);
    combined[3] = static_cast<std::byte>((blockSize >> 24U) & 0xFFU);
    std::ranges::copy(rawData, std::ranges::begin(combined) + 4);
    return combined;
}

}  // namespace

struct BamWriter::Impl
{
    BgzfWriter bgzf;
    const BamWriterConfig config;
    const bool callbacksEnabled;
    std::atomic<std::uint64_t> recordsWritten{0};
    bool closed{false};

    Impl(const std::filesystem::path& path, const SamHeader& header, const BamWriterConfig& cfg)
        : bgzf{path, cfg.BgzfConfig}, config{cfg}, callbacksEnabled{false}
    {
        const std::vector<std::byte> headerBlock{header.ToBamHeaderBlock()};
        bgzf.Write(headerBlock);
    }

    Impl(const std::filesystem::path& path, const SamHeader& header, const BamWriterConfig& cfg,
         IndexCallback cb)
        : bgzf{path, cfg.BgzfConfig}, config{cfg}, callbacksEnabled{true}
    {
        bgzf.SetCallback(std::move(cb));
        const std::vector<std::byte> headerBlock{header.ToBamHeaderBlock()};
        bgzf.Write(headerBlock);
    }
};

BamWriter::BamWriter(const std::filesystem::path& path, const SamHeader& header,
                     const BamWriterConfig& config)
    : impl_{std::make_unique<Impl>(path, header, config)}
{
}

BamWriter::BamWriter(const std::filesystem::path& path, const SamHeader& header,
                     const BamWriterConfig& config, IndexCallback callback)
    : impl_{std::make_unique<Impl>(path, header, config, std::move(callback))}
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
    std::vector<std::byte> serialized{record.SerializeToBam()};
    std::vector<std::byte> combined{BuildCombinedRecordPayload(serialized)};
    PendingCallback callback{};
    if (impl_->callbacksEnabled) {
        callback.rawData = std::move(serialized);
        callback.active = true;
    }
    impl_->bgzf.Write(std::move(combined), std::move(callback));
    impl_->recordsWritten.fetch_add(1, std::memory_order_relaxed);
}

void BamWriter::Write(std::span<const std::byte> rawData)
{
    std::vector<std::byte> combined{BuildCombinedRecordPayload(rawData)};
    PendingCallback callback{};
    if (impl_->callbacksEnabled) {
        callback.rawData =
            std::vector<std::byte>{std::ranges::begin(rawData), std::ranges::end(rawData)};
        callback.active = true;
    }
    impl_->bgzf.Write(std::move(combined), std::move(callback));
    impl_->recordsWritten.fetch_add(1, std::memory_order_relaxed);
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

WriterMetrics BamWriter::GetMetrics() const
{
    WriterMetrics m{};
    m.Bgzf = impl_->bgzf.GetMetrics();
    m.TotalRecordsWritten = impl_->recordsWritten.load(std::memory_order_relaxed);
    return m;
}

}  // namespace Samoa
}  // namespace PacBio
