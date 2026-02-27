#ifndef PBSAMOA_CORE_BGZF_HPP
#define PBSAMOA_CORE_BGZF_HPP

#include <pbsamoa/core/Metrics.hpp>
#include <pbsamoa/core/RawRecord.hpp>
#include <pbsamoa/core/SamHeader.hpp>

#include <array>
#include <compare>
#include <filesystem>
#include <functional>
#include <memory>
#include <optional>
#include <span>
#include <vector>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {

// ---------------------------------------------------------------------------
// VirtualOffset
// ---------------------------------------------------------------------------

/// \brief Type-safe virtual file offset for BGZF random access.
///
/// Encodes a compressed file offset (upper 48 bits) and an offset within
/// the decompressed block (lower 16 bits). Per the BAM spec, virtual
/// offsets support comparison but not arithmetic.
class VirtualOffset
{
public:
    constexpr VirtualOffset();
    constexpr explicit VirtualOffset(std::uint64_t rawValue);
    constexpr VirtualOffset(std::uint64_t blockOffset, std::uint16_t withinBlockOffset);

    constexpr std::uint64_t BlockOffset() const;
    constexpr std::uint16_t WithinBlockOffset() const;
    constexpr std::uint64_t Value() const;

    constexpr std::strong_ordering operator<=>(const VirtualOffset&) const = default;

private:
    std::uint64_t value_;
};

// --- inline definitions ---

constexpr VirtualOffset::VirtualOffset() : value_{0} {}

constexpr VirtualOffset::VirtualOffset(std::uint64_t rawValue) : value_{rawValue} {}

constexpr VirtualOffset::VirtualOffset(std::uint64_t blockOffset, std::uint16_t withinBlockOffset)
    : value_{(blockOffset << 16) | withinBlockOffset}
{
}

constexpr std::uint64_t VirtualOffset::BlockOffset() const { return value_ >> 16; }

constexpr std::uint16_t VirtualOffset::WithinBlockOffset() const { return value_ & 0xFFFF; }

constexpr std::uint64_t VirtualOffset::Value() const { return value_; }

// ---------------------------------------------------------------------------
// BGZF block parsing and decompression
// ---------------------------------------------------------------------------

/// \brief Parsed metadata from a BGZF block header.
struct BgzfBlockInfo
{
    std::uint32_t blockSize;             // total block size (BSIZE + 1)
    std::uint32_t compressedDataOffset;  // byte offset of CDATA within block
    std::uint32_t compressedDataSize;    // size of CDATA in bytes
};

/// \brief Parse a BGZF block header from raw bytes.
/// \param[in] data must contain at least the block header (>=18 bytes).
/// \returns block info on success, std::nullopt if not a valid BGZF header.
std::optional<BgzfBlockInfo> ParseBgzfBlockHeader(std::span<const std::byte> data);

/// \brief Decompress the CDATA portion of a BGZF block.
/// \param[in] blockData the entire BGZF block bytes
/// \param[in] info parsed header info from ParseBgzfBlockHeader
/// \param[out] output buffer to receive decompressed data (must be >= ISIZE)
/// \returns number of decompressed bytes, or nullopt on error
std::optional<std::size_t> DecompressBgzfBlock(std::span<const std::byte> blockData,
                                               BgzfBlockInfo info, std::span<std::byte> output);

/// \brief Check if 28 bytes match the BGZF EOF marker exactly.
bool IsBgzfEofMarker(std::span<const std::byte> data);

/// \brief The 28-byte EOF marker that should terminate every BGZF file.
inline constexpr std::array<std::byte, 28> BGZF_EOF_MARKER = {
    std::byte{0x1f}, std::byte{0x8b}, std::byte{0x08}, std::byte{0x04}, std::byte{0x00},
    std::byte{0x00}, std::byte{0x00}, std::byte{0x00}, std::byte{0x00}, std::byte{0xff},
    std::byte{0x06}, std::byte{0x00}, std::byte{0x42}, std::byte{0x43}, std::byte{0x02},
    std::byte{0x00}, std::byte{0x1b}, std::byte{0x00}, std::byte{0x03}, std::byte{0x00},
    std::byte{0x00}, std::byte{0x00}, std::byte{0x00}, std::byte{0x00}, std::byte{0x00},
    std::byte{0x00}, std::byte{0x00}, std::byte{0x00},
};

// ---------------------------------------------------------------------------
// BgzfReader
// ---------------------------------------------------------------------------

/// \brief Reads and decompresses BGZF-compressed files.
///
/// Two construction modes:
/// - Block-only (1-arg): synchronous ReadBlock() for index files.
/// - BAM mode (2-arg): parses BAM header eagerly and provides ReadRecord().
///   numWorkers == 0 uses synchronous record parsing.
///   numWorkers > 0 uses parallel decompression pipeline internals.
class BgzfReader
{
public:
    /// \brief Block-only mode — no header parsing, no record parsing, no threads.
    /// \throws std::runtime_error if file cannot be opened.
    explicit BgzfReader(const std::filesystem::path& path);

    /// \brief BAM mode — parses BAM header eagerly.
    /// \param[in] numWorkers 0 = synchronous record parsing, >0 = parallel pipeline.
    /// \throws std::runtime_error if file cannot be opened or header is invalid.
    BgzfReader(const std::filesystem::path& path, std::size_t numWorkers);

    ~BgzfReader();

    BgzfReader(const BgzfReader&) = delete;
    BgzfReader& operator=(const BgzfReader&) = delete;
    BgzfReader(BgzfReader&&) noexcept;
    BgzfReader& operator=(BgzfReader&&) noexcept;

    /// \brief Decompress the next BGZF block into buffer.
    /// \returns number of decompressed bytes, 0 at EOF, nullopt on error.
    std::optional<std::size_t> ReadBlock(std::span<std::byte> buffer);

    /// \brief Seek to a virtual file offset.
    void Seek(VirtualOffset offset);

    /// \brief Check if file ends with the standard 28-byte EOF marker.
    bool HasEofMarker() const;

    /// \brief Current virtual file offset (start of next unread block).
    VirtualOffset Tell() const;

    /// \brief Access the parsed BAM header.
    /// \throws std::runtime_error in block-only mode.
    const SamHeader& Header() const;

    /// \brief Read next BAM record.
    /// \throws std::runtime_error in block-only mode.
    /// \returns owning RawRecord, or nullopt at EOF.
    std::optional<RawRecord> ReadRecord();

    /// \brief Snapshot of read metrics.
    BgzfMetrics GetMetrics() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

// ---------------------------------------------------------------------------
// BgzfWriter
// ---------------------------------------------------------------------------

struct BgzfWriterConfig
{
    int CompressionLevel{6};
    std::size_t BgzfWorkers{4};
    std::size_t InputQueueCapacity{256};
    std::int32_t BlocksPerBatch{32};
};

/// \brief Metadata attached to a queued write for deferred callbacks.
struct PendingCallback
{
    std::uint16_t withinBlockOffset{0};
    std::vector<std::byte> rawData;
    std::int32_t rgId{0};
    std::int32_t zmw{0};
    bool active{false};
};

/// \brief Writes BGZF-compressed output files.
///
/// Data is accumulated into blocks of up to 64 KiB, compressed via
/// libdeflate, and written to file. An EOF marker is appended on close.
class BgzfWriter
{
public:
    using IndexCallbackFn = std::function<void(std::int64_t, std::span<const std::byte>)>;

    explicit BgzfWriter(const std::filesystem::path& path,
                        const BgzfWriterConfig& config = BgzfWriterConfig{});
    ~BgzfWriter();

    BgzfWriter(const BgzfWriter&) = delete;
    BgzfWriter& operator=(const BgzfWriter&) = delete;
    BgzfWriter(BgzfWriter&&) noexcept;
    BgzfWriter& operator=(BgzfWriter&&) noexcept;

    /// \brief Write uncompressed data. May be called multiple times.
    /// Data is buffered and compressed in 64 KiB blocks.
    void Write(std::span<const std::byte> data);

    /// \brief Write an owned payload without copying.
    void Write(std::vector<std::byte>&& data);

    /// \brief Write uncompressed data with callback metadata.
    void Write(std::span<const std::byte> data, const PendingCallback& callback);

    /// \brief Write an owned payload without copying, with callback metadata.
    void Write(std::vector<std::byte>&& data, PendingCallback callback);

    /// \brief Set callback invoked from the IO writer thread when offsets are known.
    ///
    /// Callback dispatch is asynchronous relative to Write() and may be deferred
    /// until Close() drains all queued blocks.
    void SetCallback(IndexCallbackFn callback);

    /// \brief Flush remaining buffered data and close.
    /// Called automatically by destructor.
    ///
    /// Call Close() explicitly to observe I/O errors. Errors that occur during
    /// destructor-driven close are suppressed.
    void Close();

    /// \brief Snapshot writer metrics.
    BgzfWriteMetrics GetMetrics() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_CORE_BGZF_HPP
