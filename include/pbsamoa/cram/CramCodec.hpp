#ifndef PBSAMOA_CRAM_CRAMCODEC_HPP
#define PBSAMOA_CRAM_CRAMCODEC_HPP

#include <pbsamoa/cram/CramStructs.hpp>

#include <array>
#include <memory>
#include <span>
#include <unordered_map>
#include <vector>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {

enum class CramCodecDecodeKind
{
    SCALAR,
    BYTE_ARRAY
};

// ---------------------------------------------------------------------------
// Bit reader for core data block
// ---------------------------------------------------------------------------

/// \brief Reads individual bits from a byte buffer (MSB first).
class CramBitReader
{
public:
    explicit CramBitReader(std::span<const std::byte> data);

    /// \brief Read a single bit (0 or 1).
    std::int32_t ReadBit();

    /// \brief Read numBits bits as an unsigned integer (MSB first).
    std::int32_t ReadBits(int numBits);

    /// \brief Current byte position in the buffer.
    std::size_t BytePosition() const;

private:
    std::span<const std::byte> data_;
    std::size_t bytePos_{0};
    int bitPos_{7};  // MSB first: 7..0
};

// ---------------------------------------------------------------------------
// Bit writer for core data block
// ---------------------------------------------------------------------------

/// \brief Writes individual bits to a byte buffer (MSB first).
class CramBitWriter
{
public:
    CramBitWriter();

    /// \brief Write a single bit (0 or 1).
    void WriteBit(int bit);

    /// \brief Write numBits bits from value (MSB first).
    void WriteBits(std::int32_t value, int numBits);

    /// \brief Flush remaining bits (left-shift to fill last byte).
    void Flush();

    /// \brief Get the written byte buffer (copy).
    std::vector<std::byte> Data() const&;

    /// \brief Move out the written byte buffer.
    std::vector<std::byte> Data() &&;

private:
    std::vector<std::byte> data_;
    std::byte currentByte_{};
    int bitPos_{7};
};

// ---------------------------------------------------------------------------
// External block store
// ---------------------------------------------------------------------------

/// \brief Provides byte-stream access to external data blocks by content ID.
class CramExternalBlockStore
{
public:
    /// \brief Add a block's uncompressed data by content ID.
    void AddBlock(std::int32_t contentId, std::span<const std::byte> data);

    /// \brief Read n bytes from the block with the given content ID.
    std::vector<std::byte> ReadBytes(std::int32_t contentId, std::size_t n);

    /// \brief Read n bytes as a view (no copy) from the block with the given content ID.
    std::span<const std::byte> ReadBytesView(std::int32_t contentId, std::size_t n);

    /// \brief Read bytes until a stop byte is found (zero-copy). Advances past the stop byte.
    std::span<const std::byte> ReadBytesUntilStop(std::int32_t contentId, std::byte stopByte);

    /// \brief Read a single byte from the block with the given content ID.
    std::byte ReadByte(std::int32_t contentId);

    /// \brief Read an ITF-8 integer from the block with the given content ID.
    std::int32_t ReadItf8(std::int32_t contentId);

    /// \brief Write bytes to the block with the given content ID.
    void WriteBytes(std::int32_t contentId, std::span<const std::byte> data);

    /// \brief Write a single byte to the block with the given content ID.
    void WriteByte(std::int32_t contentId, std::byte value);

    /// \brief Write an ITF-8 integer to the block with the given content ID.
    void WriteItf8(std::int32_t contentId, std::int32_t value);

    /// \brief Reserve capacity for the block with the given content ID.
    void ReserveBlock(std::int32_t contentId, std::size_t capacity);

    /// \brief Get all accumulated data for a content ID (for writing).
    std::span<const std::byte> GetBlockData(std::int32_t contentId) const;

    /// \brief Take ownership of a block's data (for writing).
    std::vector<std::byte> TakeBlockData(std::int32_t contentId);

    /// \brief Get all content IDs that have data.
    std::vector<std::int32_t> ContentIds() const;

    /// \brief Reset all read positions (for re-reading).
    void ResetPositions();

private:
    static constexpr std::size_t FAST_BLOCK_LIMIT = 256;

    struct BlockState
    {
        std::vector<std::byte> data;
        std::size_t readPos{0};
    };

    /// \brief Resolve a content ID to an existing block (throws if not found).
    BlockState* FindBlock(std::int32_t contentId);

    /// \brief Resolve a content ID to a block, creating it if necessary.
    BlockState* FindOrCreateBlock(std::int32_t contentId);

    std::array<BlockState, FAST_BLOCK_LIMIT> fastBlocks_{};
    std::array<bool, FAST_BLOCK_LIMIT> fastBlockUsed_{};
    std::unordered_map<std::int32_t, BlockState> blocks_;
};

// ---------------------------------------------------------------------------
// Codec interface
// ---------------------------------------------------------------------------

/// \brief Abstract base for CRAM data series codecs.
class CramCodec
{
public:
    CramCodec() = default;
    CramCodec(const CramCodec&) = default;
    CramCodec& operator=(const CramCodec&) = default;
    CramCodec(CramCodec&&) noexcept = default;
    CramCodec& operator=(CramCodec&&) noexcept = default;
    virtual ~CramCodec() = default;

    /// \brief Explicit payload decode mode for this codec implementation.
    virtual CramCodecDecodeKind DecodeKind() const = 0;

    /// \brief Decode a single integer value.
    virtual std::int32_t DecodeInt(CramBitReader& coreReader,
                                   CramExternalBlockStore& extBlocks) = 0;

    /// \brief Decode a single byte value.
    virtual std::byte DecodeByte(CramBitReader& coreReader, CramExternalBlockStore& extBlocks) = 0;

    /// \brief Decode a byte array.
    virtual std::vector<std::byte> DecodeByteArray(CramBitReader& coreReader,
                                                   CramExternalBlockStore& extBlocks) = 0;

    /// \brief Decode a byte array as a zero-copy view into block data.
    ///
    /// The returned span is valid until the next operation on the same block.
    /// Default implementation falls back to DecodeByteArray (with a copy).
    virtual std::span<const std::byte> DecodeByteArrayView(CramBitReader& coreReader,
                                                           CramExternalBlockStore& extBlocks);

    /// \brief Encode a single integer value.
    virtual void EncodeInt(std::int32_t value, CramBitWriter& coreWriter,
                           CramExternalBlockStore& extBlocks) = 0;

    /// \brief Encode a single byte value.
    virtual void EncodeByte(std::byte value, CramBitWriter& coreWriter,
                            CramExternalBlockStore& extBlocks) = 0;

    /// \brief Encode a byte array.
    virtual void EncodeByteArray(std::span<const std::byte> data, CramBitWriter& coreWriter,
                                 CramExternalBlockStore& extBlocks) = 0;

private:
    std::vector<std::byte> viewScratch_;
};

// ---------------------------------------------------------------------------
// Codec factory
// ---------------------------------------------------------------------------

/// \brief Create a codec instance from an encoding descriptor.
std::unique_ptr<CramCodec> CreateCodec(const CramEncodingDescriptor& desc);

// ---------------------------------------------------------------------------
// Concrete codec implementations
// ---------------------------------------------------------------------------

/// \brief NULL codec — no data is stored.
class NullCodec : public CramCodec
{
public:
    CramCodecDecodeKind DecodeKind() const override;
    std::int32_t DecodeInt(CramBitReader& coreReader, CramExternalBlockStore& extBlocks) override;
    std::byte DecodeByte(CramBitReader& coreReader, CramExternalBlockStore& extBlocks) override;
    std::vector<std::byte> DecodeByteArray(CramBitReader& coreReader,
                                           CramExternalBlockStore& extBlocks) override;
    void EncodeInt(std::int32_t value, CramBitWriter& coreWriter,
                   CramExternalBlockStore& extBlocks) override;
    void EncodeByte(std::byte value, CramBitWriter& coreWriter,
                    CramExternalBlockStore& extBlocks) override;
    void EncodeByteArray(std::span<const std::byte> data, CramBitWriter& coreWriter,
                         CramExternalBlockStore& extBlocks) override;
};

/// \brief EXTERNAL codec — data stored in an external byte-stream block.
class ExternalCodec : public CramCodec
{
public:
    explicit ExternalCodec(std::int32_t blockContentId);

    CramCodecDecodeKind DecodeKind() const override;
    std::int32_t DecodeInt(CramBitReader& coreReader, CramExternalBlockStore& extBlocks) override;
    std::byte DecodeByte(CramBitReader& coreReader, CramExternalBlockStore& extBlocks) override;
    std::vector<std::byte> DecodeByteArray(CramBitReader& coreReader,
                                           CramExternalBlockStore& extBlocks) override;
    void EncodeInt(std::int32_t value, CramBitWriter& coreWriter,
                   CramExternalBlockStore& extBlocks) override;
    void EncodeByte(std::byte value, CramBitWriter& coreWriter,
                    CramExternalBlockStore& extBlocks) override;
    void EncodeByteArray(std::span<const std::byte> data, CramBitWriter& coreWriter,
                         CramExternalBlockStore& extBlocks) override;

    /// \brief Content ID of the external block this codec reads/writes.
    std::int32_t BlockContentId() const;

private:
    std::int32_t blockContentId_;
};

/// \brief HUFFMAN codec with a single symbol (constant value).
class HuffmanCodec : public CramCodec
{
public:
    HuffmanCodec(std::vector<std::int32_t> symbols, std::vector<std::int32_t> bitLengths);

    CramCodecDecodeKind DecodeKind() const override;
    std::int32_t DecodeInt(CramBitReader& coreReader, CramExternalBlockStore& extBlocks) override;
    std::byte DecodeByte(CramBitReader& coreReader, CramExternalBlockStore& extBlocks) override;
    std::vector<std::byte> DecodeByteArray(CramBitReader& coreReader,
                                           CramExternalBlockStore& extBlocks) override;
    void EncodeInt(std::int32_t value, CramBitWriter& coreWriter,
                   CramExternalBlockStore& extBlocks) override;
    void EncodeByte(std::byte value, CramBitWriter& coreWriter,
                    CramExternalBlockStore& extBlocks) override;
    void EncodeByteArray(std::span<const std::byte> data, CramBitWriter& coreWriter,
                         CramExternalBlockStore& extBlocks) override;

private:
    struct HuffmanEntry
    {
        std::int32_t symbol;
        std::uint32_t code;
        std::int32_t bitLength;
    };

    struct HuffmanDecodeBucket
    {
        std::uint32_t minCode{0};
        std::uint32_t maxCode{0};
        std::size_t firstIndex{0};
        bool hasEntries{false};
    };

    std::vector<HuffmanEntry> entries_;
    std::vector<HuffmanDecodeBucket> decodeBuckets_;
    std::int32_t maxBitLength_{0};
    void BuildCodes();
};

/// \brief BETA codec — fixed-width binary encoding.
class BetaCodec : public CramCodec
{
public:
    BetaCodec(std::int32_t offset, std::int32_t numBits);

    CramCodecDecodeKind DecodeKind() const override;
    std::int32_t DecodeInt(CramBitReader& coreReader, CramExternalBlockStore& extBlocks) override;
    std::byte DecodeByte(CramBitReader& coreReader, CramExternalBlockStore& extBlocks) override;
    std::vector<std::byte> DecodeByteArray(CramBitReader& coreReader,
                                           CramExternalBlockStore& extBlocks) override;
    void EncodeInt(std::int32_t value, CramBitWriter& coreWriter,
                   CramExternalBlockStore& extBlocks) override;
    void EncodeByte(std::byte value, CramBitWriter& coreWriter,
                    CramExternalBlockStore& extBlocks) override;
    void EncodeByteArray(std::span<const std::byte> data, CramBitWriter& coreWriter,
                         CramExternalBlockStore& extBlocks) override;

private:
    std::int32_t offset_;
    std::int32_t numBits_;
};

/// \brief GAMMA codec — Elias gamma encoding.
class GammaCodec : public CramCodec
{
public:
    explicit GammaCodec(std::int32_t offset);

    CramCodecDecodeKind DecodeKind() const override;
    std::int32_t DecodeInt(CramBitReader& coreReader, CramExternalBlockStore& extBlocks) override;
    std::byte DecodeByte(CramBitReader& coreReader, CramExternalBlockStore& extBlocks) override;
    std::vector<std::byte> DecodeByteArray(CramBitReader& coreReader,
                                           CramExternalBlockStore& extBlocks) override;
    void EncodeInt(std::int32_t value, CramBitWriter& coreWriter,
                   CramExternalBlockStore& extBlocks) override;
    void EncodeByte(std::byte value, CramBitWriter& coreWriter,
                    CramExternalBlockStore& extBlocks) override;
    void EncodeByteArray(std::span<const std::byte> data, CramBitWriter& coreWriter,
                         CramExternalBlockStore& extBlocks) override;

private:
    std::int32_t offset_;
};

/// \brief SUBEXP codec — subexponential encoding.
class SubexpCodec : public CramCodec
{
public:
    SubexpCodec(std::int32_t offset, std::int32_t k);

    CramCodecDecodeKind DecodeKind() const override;
    std::int32_t DecodeInt(CramBitReader& coreReader, CramExternalBlockStore& extBlocks) override;
    std::byte DecodeByte(CramBitReader& coreReader, CramExternalBlockStore& extBlocks) override;
    std::vector<std::byte> DecodeByteArray(CramBitReader& coreReader,
                                           CramExternalBlockStore& extBlocks) override;
    void EncodeInt(std::int32_t value, CramBitWriter& coreWriter,
                   CramExternalBlockStore& extBlocks) override;
    void EncodeByte(std::byte value, CramBitWriter& coreWriter,
                    CramExternalBlockStore& extBlocks) override;
    void EncodeByteArray(std::span<const std::byte> data, CramBitWriter& coreWriter,
                         CramExternalBlockStore& extBlocks) override;

private:
    std::int32_t offset_;
    std::int32_t k_;
};

/// \brief BYTE_ARRAY_LEN codec — byte array with explicit length.
class ByteArrayLenCodec : public CramCodec
{
public:
    ByteArrayLenCodec(std::unique_ptr<CramCodec> lenCodec, std::unique_ptr<CramCodec> dataCodec);

    CramCodecDecodeKind DecodeKind() const override;
    std::int32_t DecodeInt(CramBitReader& coreReader, CramExternalBlockStore& extBlocks) override;
    std::byte DecodeByte(CramBitReader& coreReader, CramExternalBlockStore& extBlocks) override;
    std::vector<std::byte> DecodeByteArray(CramBitReader& coreReader,
                                           CramExternalBlockStore& extBlocks) override;
    std::span<const std::byte> DecodeByteArrayView(CramBitReader& coreReader,
                                                   CramExternalBlockStore& extBlocks) override;
    void EncodeInt(std::int32_t value, CramBitWriter& coreWriter,
                   CramExternalBlockStore& extBlocks) override;
    void EncodeByte(std::byte value, CramBitWriter& coreWriter,
                    CramExternalBlockStore& extBlocks) override;
    void EncodeByteArray(std::span<const std::byte> data, CramBitWriter& coreWriter,
                         CramExternalBlockStore& extBlocks) override;

private:
    std::unique_ptr<CramCodec> lenCodec_;
    std::unique_ptr<CramCodec> dataCodec_;
    std::vector<std::byte> viewScratch2_;
};

/// \brief BYTE_ARRAY_STOP codec — byte array terminated by stop byte.
class ByteArrayStopCodec : public CramCodec
{
public:
    ByteArrayStopCodec(std::byte stopByte, std::int32_t blockContentId);

    CramCodecDecodeKind DecodeKind() const override;
    std::int32_t DecodeInt(CramBitReader& coreReader, CramExternalBlockStore& extBlocks) override;
    std::byte DecodeByte(CramBitReader& coreReader, CramExternalBlockStore& extBlocks) override;
    std::vector<std::byte> DecodeByteArray(CramBitReader& coreReader,
                                           CramExternalBlockStore& extBlocks) override;
    std::span<const std::byte> DecodeByteArrayView(CramBitReader& coreReader,
                                                   CramExternalBlockStore& extBlocks) override;
    void EncodeInt(std::int32_t value, CramBitWriter& coreWriter,
                   CramExternalBlockStore& extBlocks) override;
    void EncodeByte(std::byte value, CramBitWriter& coreWriter,
                    CramExternalBlockStore& extBlocks) override;
    void EncodeByteArray(std::span<const std::byte> data, CramBitWriter& coreWriter,
                         CramExternalBlockStore& extBlocks) override;

private:
    std::byte stopByte_;
    std::int32_t blockContentId_;
};

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_CRAM_CRAMCODEC_HPP
