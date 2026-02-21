#ifndef PBSAMOA_IO_BAMWRITER_HPP
#define PBSAMOA_IO_BAMWRITER_HPP

#include <pbsamoa/core/BamRecord.hpp>
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

/// \brief Callback invoked after each record write with virtual offset and raw data.
using IndexCallback =
    std::function<void(std::int64_t virtualOffset, std::span<const std::byte> rawData)>;

/// \brief Writes BAM files via BGZF compression.
///
/// Accepts owned BamRecord (serializes to BAM binary), raw byte spans
/// (zero-copy passthrough), or RawRecord (delegates to raw path).
class BamWriter
{
public:
    /// \param[in] path output BAM file path
    /// \param[in] header SAM header to write
    /// \param[in] compressionLevel libdeflate level 1-12 (default 6)
    BamWriter(const std::filesystem::path& path, const SamHeader& header, int compressionLevel = 6);

    /// \param[in] path output BAM file path
    /// \param[in] header SAM header to write
    /// \param[in] compressionLevel libdeflate level 1-12
    /// \param[in] callback invoked after each record write with virtual offset and raw data
    BamWriter(const std::filesystem::path& path, const SamHeader& header, int compressionLevel,
              IndexCallback callback);
    ~BamWriter();

    BamWriter(const BamWriter&) = delete;
    BamWriter& operator=(const BamWriter&) = delete;
    BamWriter(BamWriter&&) noexcept;
    BamWriter& operator=(BamWriter&&) noexcept;

    /// \brief Write an owned record (serializes to BAM binary).
    void Write(const BamRecord& record);

    /// \brief Write raw record bytes (zero-copy passthrough).
    void Write(std::span<const std::byte> rawData);

    /// \brief Write an owning byte view (delegates to raw path).
    void Write(const RawRecord& byteView);

    /// \brief Write a batch of records (raw zero-copy path).
    void WriteBatch(const RawRecordBatch& batch);

    /// \brief Flush and close. Called automatically by destructor.
    void Close();

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_IO_BAMWRITER_HPP
