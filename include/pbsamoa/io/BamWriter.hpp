#ifndef PBSAMOA_IO_BAMWRITER_HPP
#define PBSAMOA_IO_BAMWRITER_HPP

#include <pbsamoa/core/BamRecord.hpp>
#include <pbsamoa/core/Bgzf.hpp>
#include <pbsamoa/core/RawRecord.hpp>
#include <pbsamoa/core/SamHeader.hpp>

#include <filesystem>
#include <functional>
#include <memory>
#include <span>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {

/// \brief Callback invoked with virtual offset and raw data for each record.
///
/// Callbacks are dispatched asynchronously from BgzfWriter's IO thread once
/// final block offsets are known (possibly during Close()).
using IndexCallback =
    std::function<void(std::int64_t virtualOffset, std::span<const std::byte> rawData)>;

struct BamWriterConfig
{
    BgzfWriterConfig BgzfConfig{};
    bool UseTempFile{false};
};

/// \brief Writes BAM files via BGZF compression.
///
/// Accepts owned BamRecord (serializes to BAM binary), raw byte spans,
/// RawRecordView, or RawRecord (all delegate to the raw path).
class BamWriter
{
public:
    BamWriter(const std::filesystem::path& path, const SamHeader& header,
              const BamWriterConfig& config = BamWriterConfig{});

    BamWriter(const std::filesystem::path& path, const SamHeader& header,
              const BamWriterConfig& config, IndexCallback callback);
    ~BamWriter();

    BamWriter(const BamWriter&) = delete;
    BamWriter& operator=(const BamWriter&) = delete;
    BamWriter(BamWriter&&) noexcept;
    BamWriter& operator=(BamWriter&&) noexcept;

    /// \brief Write an owned record (serializes to BAM binary).
    void Write(const BamRecord& record);

    /// \brief Write raw record bytes (zero-copy passthrough).
    void Write(std::span<const std::byte> rawData);

    /// \brief Write a non-owning raw record view (delegates to raw path).
    void Write(const RawRecordView& byteView);

    /// \brief Write an owning raw record (delegates to raw path).
    void Write(const RawRecord& byteView);

    /// \brief Write a batch of records (raw zero-copy path).
    void WriteBatch(const RawRecordBatch& batch);

    /// \brief Flush and close. Called automatically by destructor.
    void Close();

    /// \brief Snapshot writer metrics.
    WriterMetrics GetMetrics() const;

    /// \brief Virtual offset just past the last record written (for on-the-fly indexing).
    ///
    /// Pairs with the per-record index callback: the callback supplies each record's
    /// begin offset, and this supplies the end offset of the final record. Meaningful
    /// only after Close().
    VirtualOffset EndVirtualOffset() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

/// \brief Append one record (block_size field and record bytes) to \p encoder.
///
/// Use it with BgzfBlockWriter to serialize and compress records on many threads:
/// each thread fills its own encoder, and the writer only appends the blocks.
/// Write the BAM header block (SamHeader::ToBamHeaderBlock()) first.
/// \throws std::invalid_argument if the record exceeds the BAM block_size limit
void EncodeBamRecord(const BamRecord& record, BgzfBlockEncoder& encoder);

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_IO_BAMWRITER_HPP
