#ifndef PBSAMOA_CRAM_CRAMSTRUCTS_HPP
#define PBSAMOA_CRAM_CRAMSTRUCTS_HPP

#include <array>
#include <span>
#include <vector>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {

namespace detail {

inline constexpr std::uint16_t MakeCramDataSeriesCode(const char first, const char second)
{
    return (static_cast<std::uint16_t>(static_cast<std::uint8_t>(first)) << 8) |
           static_cast<std::uint16_t>(static_cast<std::uint8_t>(second));
}

}  // namespace detail

// ---------------------------------------------------------------------------
// ITF-8 / LTF-8 variable-length integer encoding
// ---------------------------------------------------------------------------

/// \brief Read an ITF-8 encoded integer from a byte span.
/// \param[in] data input bytes
/// \param[out] bytesRead number of bytes consumed
/// \returns decoded 32-bit value
std::int32_t ReadItf8(std::span<const std::byte> data, std::size_t& bytesRead);

/// \brief Read an LTF-8 encoded long from a byte span.
/// \param[in] data input bytes
/// \param[out] bytesRead number of bytes consumed
/// \returns decoded 64-bit value
std::int64_t ReadLtf8(std::span<const std::byte> data, std::size_t& bytesRead);

/// \brief Write an ITF-8 encoded integer to a byte vector.
/// \returns number of bytes written
std::size_t WriteItf8(std::vector<std::byte>& out, std::int32_t value);

/// \brief Write an LTF-8 encoded long to a byte vector.
/// \returns number of bytes written
std::size_t WriteLtf8(std::vector<std::byte>& out, std::int64_t value);

// ---------------------------------------------------------------------------
// File Definition
// ---------------------------------------------------------------------------

/// \brief CRAM magic bytes: "CRAM"
inline constexpr std::array<std::byte, 4> CRAM_MAGIC = {std::byte{0x43}, std::byte{0x52},
                                                        std::byte{0x41}, std::byte{0x4d}};

struct CramFileDefinition
{
    std::uint8_t MajorVersion{3};
    std::uint8_t MinorVersion{0};
    std::array<std::byte, 20> FileId{};
};

// ---------------------------------------------------------------------------
// Container Header
// ---------------------------------------------------------------------------

struct CramContainerHeader
{
    std::int32_t Length{0};
    std::int32_t RefSeqId{0};
    std::int32_t StartPos{0};
    std::int32_t AlignmentSpan{0};
    std::int32_t NumRecords{0};
    std::int64_t RecordCounter{0};
    std::int64_t Bases{0};
    std::int32_t NumBlocks{0};
    std::vector<std::int32_t> Landmarks;
    std::uint32_t Crc32{0};
};

// ---------------------------------------------------------------------------
// Block
// ---------------------------------------------------------------------------

/// \brief Block compression methods.
enum class CramBlockMethod : std::uint8_t
{
    RAW = 0,
    GZIP = 1,
    BZIP2 = 2,
    LZMA = 3,
    RANS4X8 = 4,
    RANS4X16 = 5,
    ADAPTIVE_ARITH = 6,
    FQZCOMP = 7,
    NAME_TOKENISER = 8,
};

/// \brief Block content types.
enum class CramBlockContentType : std::uint8_t
{
    FILE_HEADER = 0,
    COMPRESSION_HEADER = 1,
    SLICE_HEADER = 2,
    RESERVED = 3,
    EXTERNAL_DATA = 4,
    CORE_DATA = 5,
};

struct CramBlock
{
    CramBlockMethod Method{CramBlockMethod::RAW};
    CramBlockContentType ContentType{CramBlockContentType::FILE_HEADER};
    std::int32_t ContentId{0};
    std::int32_t CompressedSize{0};
    std::int32_t RawSize{0};
    std::vector<std::byte> Data;  // uncompressed data after reading; raw data for writing
    std::uint32_t Crc32{0};
};

// ---------------------------------------------------------------------------
// Slice Header
// ---------------------------------------------------------------------------

struct CramSliceHeader
{
    std::int32_t RefSeqId{0};
    std::int32_t AlignmentStart{0};
    std::int32_t AlignmentSpan{0};
    std::int32_t NumRecords{0};
    std::int64_t RecordCounter{0};
    std::int32_t NumBlocks{0};
    std::vector<std::int32_t> BlockContentIds;
    std::int32_t EmbeddedRefBlockId{-1};
    std::array<std::byte, 16> RefMd5{};
    std::vector<std::byte> OptionalTags;  // BAM-encoded optional tags
};

// ---------------------------------------------------------------------------
// Preservation Map
// ---------------------------------------------------------------------------

struct CramPreservationMap
{
    bool ReadNamesIncluded{true};
    bool ApDelta{true};
    bool ReferenceRequired{true};
    // Default CRAM substitution matrix equivalent to "CGTNAGTNACTNACGNACGT".
    std::array<std::byte, 5> SubstitutionMatrix{std::byte{0x1B}, std::byte{0x1B}, std::byte{0x1B},
                                                std::byte{0x1B}, std::byte{0x1B}};
    std::vector<std::byte> TagIdsDictionary;
};

// ---------------------------------------------------------------------------
// Encoding descriptor
// ---------------------------------------------------------------------------

/// \brief Codec IDs used in the data series encoding map.
enum class CramCodecId : std::int32_t
{
    NONE = 0,
    EXTERNAL = 1,
    HUFFMAN = 3,
    BYTE_ARRAY_LEN = 4,
    BYTE_ARRAY_STOP = 5,
    BETA = 6,
    SUBEXP = 7,
    GAMMA = 9,
};

/// \brief A serialized encoding descriptor: codec ID + raw parameter bytes.
struct CramEncodingDescriptor
{
    CramCodecId CodecId{CramCodecId::NONE};
    std::vector<std::byte> Parameters;
};

// ---------------------------------------------------------------------------
// Data Series Keys (2-byte codes)
// ---------------------------------------------------------------------------

/// \brief Known CRAM data series keys.
enum class CramDataSeries : std::uint16_t
{
    BF = detail::MakeCramDataSeriesCode('B', 'F'),
    CF = detail::MakeCramDataSeriesCode('C', 'F'),
    RI = detail::MakeCramDataSeriesCode('R', 'I'),
    RL = detail::MakeCramDataSeriesCode('R', 'L'),
    AP = detail::MakeCramDataSeriesCode('A', 'P'),
    RG = detail::MakeCramDataSeriesCode('R', 'G'),
    RN = detail::MakeCramDataSeriesCode('R', 'N'),
    MF = detail::MakeCramDataSeriesCode('M', 'F'),
    NS = detail::MakeCramDataSeriesCode('N', 'S'),
    NP = detail::MakeCramDataSeriesCode('N', 'P'),
    TS = detail::MakeCramDataSeriesCode('T', 'S'),
    NF = detail::MakeCramDataSeriesCode('N', 'F'),
    TL = detail::MakeCramDataSeriesCode('T', 'L'),
    FN = detail::MakeCramDataSeriesCode('F', 'N'),
    FC = detail::MakeCramDataSeriesCode('F', 'C'),
    FP = detail::MakeCramDataSeriesCode('F', 'P'),
    DL = detail::MakeCramDataSeriesCode('D', 'L'),
    BB = detail::MakeCramDataSeriesCode('B', 'B'),
    QQ = detail::MakeCramDataSeriesCode('Q', 'Q'),
    BS = detail::MakeCramDataSeriesCode('B', 'S'),
    IN = detail::MakeCramDataSeriesCode('I', 'N'),
    RS = detail::MakeCramDataSeriesCode('R', 'S'),
    PD = detail::MakeCramDataSeriesCode('P', 'D'),
    HC = detail::MakeCramDataSeriesCode('H', 'C'),
    SC = detail::MakeCramDataSeriesCode('S', 'C'),
    MQ = detail::MakeCramDataSeriesCode('M', 'Q'),
    BA = detail::MakeCramDataSeriesCode('B', 'A'),
    QS = detail::MakeCramDataSeriesCode('Q', 'S'),
};

// ---------------------------------------------------------------------------
// Compression Header
// ---------------------------------------------------------------------------

struct CramCompressionHeader
{
    CramPreservationMap PreservationMap;
    std::vector<std::pair<CramDataSeries, CramEncodingDescriptor>> DataSeriesEncodings;
    std::vector<std::pair<std::int32_t, CramEncodingDescriptor>> TagEncodings;
};

// ---------------------------------------------------------------------------
// Composed Types
// ---------------------------------------------------------------------------

/// \brief A fully-formed CRAM slice: parsed slice header + data blocks.
struct CramSlice
{
    CramSliceHeader Header;
    CramBlock CoreBlock;
    std::vector<CramBlock> ExternalBlocks;
};

/// \brief A fully-formed CRAM container: header + compression header + slices.
struct CramContainer
{
    CramContainerHeader Header;
    CramCompressionHeader CompressionHeader;
    std::vector<CramSlice> Slices;
};

// ---------------------------------------------------------------------------
// EOF marker
// ---------------------------------------------------------------------------

/// \brief The 38-byte CRAM EOF container marker.
inline constexpr std::array<std::byte, 38> CRAM_EOF_MARKER = {
    std::byte{0x0f}, std::byte{0x00}, std::byte{0x00}, std::byte{0x00}, std::byte{0xff},
    std::byte{0xff}, std::byte{0xff}, std::byte{0xff}, std::byte{0x0f}, std::byte{0xe0},
    std::byte{0x45}, std::byte{0x4f}, std::byte{0x46}, std::byte{0x00}, std::byte{0x00},
    std::byte{0x00}, std::byte{0x00}, std::byte{0x01}, std::byte{0x00}, std::byte{0x05},
    std::byte{0xbd}, std::byte{0xd9}, std::byte{0x4f}, std::byte{0x00}, std::byte{0x01},
    std::byte{0x00}, std::byte{0x06}, std::byte{0x06}, std::byte{0x01}, std::byte{0x00},
    std::byte{0x01}, std::byte{0x00}, std::byte{0x01}, std::byte{0x00}, std::byte{0xee},
    std::byte{0x63}, std::byte{0x01}, std::byte{0x4b},
};

// ---------------------------------------------------------------------------
// Parsing functions
// ---------------------------------------------------------------------------

/// \brief Read a CRAM file definition from raw bytes.
CramFileDefinition ParseFileDefinition(std::span<const std::byte> data);

/// \brief Serialize a CRAM file definition to bytes.
std::vector<std::byte> SerializeFileDefinition(const CramFileDefinition& def);

/// \brief Read a container header from a byte span.
/// \param[out] bytesRead number of bytes consumed
CramContainerHeader ParseContainerHeader(std::span<const std::byte> data, std::size_t& bytesRead);

/// \brief Serialize a container header to bytes.
std::vector<std::byte> SerializeContainerHeader(const CramContainerHeader& header);

/// \brief Read a block (header + data + CRC32) from a byte span.
/// \param[out] bytesRead number of bytes consumed
CramBlock ParseBlock(std::span<const std::byte> data, std::size_t& bytesRead);

/// \brief Serialize a block to bytes (header + compressed data + CRC32).
std::vector<std::byte> SerializeBlock(const CramBlock& block);

/// \brief Parse a slice header from uncompressed block data.
CramSliceHeader ParseSliceHeader(std::span<const std::byte> data);

/// \brief Serialize a slice header to bytes.
std::vector<std::byte> SerializeSliceHeader(const CramSliceHeader& header);

/// \brief Parse a compression header from uncompressed block data.
CramCompressionHeader ParseCompressionHeader(std::span<const std::byte> data);

/// \brief Serialize a compression header to bytes.
std::vector<std::byte> SerializeCompressionHeader(const CramCompressionHeader& header);

/// \brief Serialize a complete slice (slice header block + data blocks).
std::vector<std::byte> SerializeSlice(const CramSlice& slice);

/// \brief Compute serialized byte size of a slice without allocating output
/// bytes.
std::size_t SerializedSliceSize(const CramSlice& slice);

/// \brief Serialize a complete container (header + payload blocks).
std::vector<std::byte> SerializeContainer(const CramContainer& container);

/// \brief Parse a full container payload (bytes after the container header).
///
/// \param header Parsed container header for this payload.
/// \param payload Raw payload bytes of length \c header.Length.
/// \returns Parsed container with decompressed compression/slice headers and
///          still-compressed data blocks.
CramContainer ParseContainer(const CramContainerHeader& header, std::span<const std::byte> payload);

/// \brief Check if bytes match the CRAM EOF marker.
bool IsCramEofMarker(std::span<const std::byte> data);

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_CRAM_CRAMSTRUCTS_HPP
