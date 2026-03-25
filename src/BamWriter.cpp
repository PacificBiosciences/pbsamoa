#include <pbsamoa/io/BamWriter.hpp>

#include <algorithm>
#include <atomic>
#include <format>
#include <limits>
#include <stdexcept>
#include <vector>

#include <cstdint>

namespace PacBio {
namespace Samoa {

namespace {

constexpr std::size_t BAM_BLOCK_SIZE_PREFIX_BYTES{4};

std::vector<std::byte> BuildCombinedRecordPayload(std::span<const std::byte> rawData)
{
    if (std::size(rawData) > static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max())) {
        throw std::invalid_argument{
            std::format("BamWriter::Write: raw record size ({}) exceeds BAM block_size limit",
                        std::size(rawData))};
    }

    const std::uint32_t blockSize{static_cast<std::uint32_t>(std::size(rawData))};

    std::vector<std::byte> combined{};
    combined.resize(BAM_BLOCK_SIZE_PREFIX_BYTES + std::size(rawData));
    combined[0] = static_cast<std::byte>(blockSize & 0xFFU);
    combined[1] = static_cast<std::byte>((blockSize >> 8U) & 0xFFU);
    combined[2] = static_cast<std::byte>((blockSize >> 16U) & 0xFFU);
    combined[3] = static_cast<std::byte>((blockSize >> 24U) & 0xFFU);
    std::ranges::copy(rawData, std::ranges::begin(combined) + BAM_BLOCK_SIZE_PREFIX_BYTES);
    return combined;
}

void EmitRawRecord(BgzfWriter& bgzf, std::atomic<std::uint64_t>& recordsWritten,
                   std::vector<std::byte> combined, PendingCallback callback)
{
    bgzf.Write(std::move(combined), std::move(callback));
    recordsWritten.fetch_add(1, std::memory_order_relaxed);
}

void WriteRawRecord(BgzfWriter& bgzf, std::atomic<std::uint64_t>& recordsWritten,
                    std::span<const std::byte> rawData, PendingCallback callback)
{
    EmitRawRecord(bgzf, recordsWritten, BuildCombinedRecordPayload(rawData), std::move(callback));
}

void WriteRawRecord(BgzfWriter& bgzf, std::atomic<std::uint64_t>& recordsWritten,
                    std::vector<std::byte>&& rawData)
{
    PendingCallback callback{};
    callback.rawData = std::move(rawData);
    callback.active = true;
    EmitRawRecord(bgzf, recordsWritten, BuildCombinedRecordPayload(callback.rawData),
                  std::move(callback));
}

void WriteHeaderBlock(BgzfWriter& bgzf, const SamHeader& header)
{
    const std::vector<std::byte> headerBlock{header.ToBamHeaderBlock()};
    bgzf.Write(headerBlock);
}

}  // namespace

struct BamWriter::Impl
{
    BgzfWriter bgzf;
    bool callbacksEnabled{false};
    std::atomic<std::uint64_t> recordsWritten{0};
    bool closed{false};

    static BgzfWriterConfig MergeBgzfConfig(const BamWriterConfig& cfg)
    {
        BgzfWriterConfig bgzfConfig{cfg.BgzfConfig};
        bgzfConfig.UseTempFile = bgzfConfig.UseTempFile || cfg.UseTempFile;
        return bgzfConfig;
    }

    Impl(const std::filesystem::path& path, const SamHeader& header, const BamWriterConfig& cfg)
        : bgzf{path, MergeBgzfConfig(cfg)}, callbacksEnabled{false}
    {
        WriteHeaderBlock(bgzf, header);
    }

    Impl(const std::filesystem::path& path, const SamHeader& header, const BamWriterConfig& cfg,
         IndexCallback cb)
        : bgzf{path, MergeBgzfConfig(cfg)}, callbacksEnabled{true}
    {
        bgzf.SetCallback(std::move(cb));
        WriteHeaderBlock(bgzf, header);
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
    if (impl_ && !impl_->closed) {
        Close();
    }
}

BamWriter::BamWriter(BamWriter&&) noexcept = default;
BamWriter& BamWriter::operator=(BamWriter&&) noexcept = default;

void BamWriter::Write(const BamRecord& record)
{
    std::vector<std::byte> serialized{record.SerializeToBam()};
    if (impl_->callbacksEnabled) {
        WriteRawRecord(impl_->bgzf, impl_->recordsWritten, std::move(serialized));
        return;
    }

    WriteRawRecord(impl_->bgzf, impl_->recordsWritten, std::span<const std::byte>{serialized},
                   PendingCallback{});
}

void BamWriter::Write(std::span<const std::byte> rawData)
{
    if (impl_->callbacksEnabled) {
        std::vector<std::byte> rawBytes{std::ranges::begin(rawData), std::ranges::end(rawData)};
        WriteRawRecord(impl_->bgzf, impl_->recordsWritten, std::move(rawBytes));
        return;
    }

    WriteRawRecord(impl_->bgzf, impl_->recordsWritten, rawData, PendingCallback{});
}

void BamWriter::Write(const RawRecord& byteView) { Write(byteView.RawData()); }

void BamWriter::WriteBatch(const RawRecordBatch& batch)
{
    const std::size_t recordCount{batch.RecordCount()};
    for (std::size_t i{0}; i < recordCount; ++i) {
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
