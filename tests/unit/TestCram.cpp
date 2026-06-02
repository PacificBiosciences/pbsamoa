#include "TestData.hpp"
#include "TestTempDir.hpp"

#include "../../src/BinaryUtils.hpp"
#include "../../src/CramInternal.hpp"

#include <pbsamoa/cram/CramCodec.hpp>
#include <pbsamoa/cram/CramCompression.hpp>
#include <pbsamoa/cram/CramStructs.hpp>
#include <pbsamoa/index/CraiIndex.hpp>
#include <pbsamoa/io/CramReader.hpp>
#include <pbsamoa/io/CramWriter.hpp>

#include <pbsamoa/core/BamRecord.hpp>
#include <pbsamoa/core/CigarOp.hpp>

#include <gtest/gtest.h>
#include <libdeflate.h>

#include <array>
#include <filesystem>
#include <format>
#include <fstream>
#include <memory>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {

namespace {

std::filesystem::path DefaultCraiPathFor(const std::filesystem::path& cramPath)
{
    return std::filesystem::path{cramPath.string() + ".crai"};
}

std::size_t FqzcompRoundTripRecordCount(const std::int32_t recordsPerSlice)
{
    if (recordsPerSlice == 1000) {
        return 1000;
    }
    return 37;
}

std::int32_t QueryRecordEnd(const BamRecord& record)
{
    const std::int32_t refEnd = record.ReferenceEnd();
    if (refEnd > record.Pos()) {
        return refEnd;
    }
    return record.Pos() + 1;
}

std::vector<std::byte> ReadFileBytes(const std::filesystem::path& path)
{
    std::ifstream in{path, std::ios::binary};
    if (!in.is_open()) {
        throw std::runtime_error("failed to open file: " + path.string());
    }

    in.seekg(0, std::ios::end);
    const std::streamoff fileSize = in.tellg();
    if (fileSize < 0) {
        throw std::runtime_error("failed to read file size: " + path.string());
    }
    in.seekg(0, std::ios::beg);

    std::vector<std::byte> bytes(static_cast<std::size_t>(fileSize));
    if (fileSize > 0) {
        in.read(reinterpret_cast<char*>(bytes.data()), fileSize);
        if (in.gcount() != fileSize) {
            throw std::runtime_error("failed to read file: " + path.string());
        }
    }
    return bytes;
}

struct CramSliceOffsetInfo
{
    std::int64_t ContainerOffset{};
    std::int64_t SliceOffset{};
    std::int64_t SliceSize{};
};

CramSliceOffsetInfo FirstDataSliceOffsetInfo(const std::filesystem::path& cramPath)
{
    constexpr std::size_t CRAM_FILE_DEFINITION_SIZE = 26;
    const std::vector<std::byte> bytes = ReadFileBytes(cramPath);
    if (std::size(bytes) <= CRAM_FILE_DEFINITION_SIZE) {
        throw std::runtime_error("CRAM too small to contain a data container");
    }

    std::size_t offset = CRAM_FILE_DEFINITION_SIZE;
    std::size_t headerContainerBytes = 0;
    const auto headerContainer = ParseContainerHeader(
        std::span<const std::byte>{bytes.data() + offset, std::size(bytes) - offset},
        headerContainerBytes);
    offset += headerContainerBytes + static_cast<std::size_t>(headerContainer.Length);

    std::size_t dataContainerHeaderBytes = 0;
    const auto dataContainer = ParseContainerHeader(
        std::span<const std::byte>{bytes.data() + offset, std::size(bytes) - offset},
        dataContainerHeaderBytes);
    if (std::empty(dataContainer.Landmarks)) {
        throw std::runtime_error("CRAM data container missing landmarks");
    }

    return CramSliceOffsetInfo{
        .ContainerOffset = static_cast<std::int64_t>(offset),
        .SliceOffset = dataContainer.Landmarks.front(),
        .SliceSize =
            static_cast<std::int64_t>(dataContainer.Length) - dataContainer.Landmarks.front(),
    };
}

CramSlice MakeRawExternalSlice(const std::int32_t recordCounter, const std::int32_t contentId)
{
    CramSlice slice;
    slice.Header.RefSeqId = -1;
    slice.Header.AlignmentStart = 1;
    slice.Header.AlignmentSpan = 0;
    slice.Header.NumRecords = 1;
    slice.Header.RecordCounter = recordCounter;
    slice.Header.NumBlocks = 2;
    slice.Header.BlockContentIds = {contentId};
    slice.Header.EmbeddedRefBlockId = -1;

    slice.CoreBlock.Method = CramBlockMethod::RAW;
    slice.CoreBlock.ContentType = CramBlockContentType::CORE_DATA;
    slice.CoreBlock.ContentId = 0;
    slice.CoreBlock.RawSize = 0;
    slice.CoreBlock.CompressedSize = 0;
    slice.CoreBlock.Data = {};

    CramBlock ext;
    ext.Method = CramBlockMethod::RAW;
    ext.ContentType = CramBlockContentType::EXTERNAL_DATA;
    ext.ContentId = contentId;
    ext.Data = {std::byte{0x01}, std::byte{0x02}};
    ext.RawSize = 2;
    ext.CompressedSize = 2;
    slice.ExternalBlocks.push_back(std::move(ext));
    return slice;
}

void AppendSerializedRecord(std::vector<std::byte>& buffer,
                            std::vector<RawRecordBatch::RecordExtent>& extents,
                            const BamRecord& record)
{
    const std::vector<std::byte> bytes = record.SerializeToBam();
    const auto offset = static_cast<std::uint32_t>(std::size(buffer));
    buffer.insert(std::end(buffer), std::begin(bytes), std::end(bytes));
    extents.push_back(
        RawRecordBatch::RecordExtent{offset, static_cast<std::uint32_t>(std::size(bytes))});
}

void ExpectDecodedIntTag(const BamRecord& record, const TagKey key, const std::int64_t expected)
{
    const auto* tag = record.Tags().Get(key);
    ASSERT_TRUE(tag);
    const auto* value = std::get_if<std::int64_t>(tag);
    ASSERT_TRUE(value);
    EXPECT_EQ(*value, expected);
}

std::vector<std::byte> EmptyCompressionResult(std::span<const std::byte> /*data*/) { return {}; }

std::vector<std::byte> EmptyDecompressionResult(std::span<const std::byte> /*data*/,
                                                std::size_t /*rawSize*/)
{
    return {};
}

std::vector<std::byte> XorCompress(std::span<const std::byte> data)
{
    std::vector<std::byte> result(data.begin(), data.end());
    for (auto& byte : result) {
        byte ^= std::byte{0xFF};
    }
    return result;
}

std::vector<std::byte> XorDecompress(std::span<const std::byte> data, std::size_t /*rawSize*/)
{
    return XorCompress(data);
}

void WriteSingleRecord(const std::filesystem::path& path, const SamHeader& header,
                       const CramWriterConfig& config, const BamRecord& record)
{
    CramWriter writer{path, header, config};
    writer.Write(record);
    writer.Close();
}

}  // namespace

// ===========================================================================
// ITF-8 / LTF-8 Tests
// ===========================================================================

TEST(CramInternal, TagContentIdMatchesPackedTagBytes)
{
    constexpr TagKey key{'N', 'M'};
    constexpr char type = 'i';
    constexpr std::int32_t expected =
        (static_cast<std::int32_t>(static_cast<std::uint8_t>('N')) << 16) |
        (static_cast<std::int32_t>(static_cast<std::uint8_t>('M')) << 8) |
        static_cast<std::int32_t>(static_cast<std::uint8_t>('i'));

    EXPECT_EQ(TagContentId('N', 'M', type), expected);
    EXPECT_EQ(TagContentId(key, type), expected);
}

TEST(CramInternal, TagTripleOrderingIsStable)
{
    std::vector<TagTriple> values{
        TagTriple{'R', 'G', 'Z'},
        TagTriple{'A', 'A', 'i'},
        TagTriple{'A', 'A', 'Z'},
        TagTriple{'A', 'B', 'i'},
    };

    std::ranges::sort(values);

    ASSERT_EQ(values.size(), 4u);
    EXPECT_EQ(values[0], (TagTriple{'A', 'A', 'Z'}));
    EXPECT_EQ(values[1], (TagTriple{'A', 'A', 'i'}));
    EXPECT_EQ(values[2], (TagTriple{'A', 'B', 'i'}));
    EXPECT_EQ(values[3], (TagTriple{'R', 'G', 'Z'}));
}

TEST(CramInternal, ReadI32LEReadsLittleEndianValues)
{
    const std::array<std::byte, 4> bytes{
        std::byte{0x78},
        std::byte{0x56},
        std::byte{0x34},
        std::byte{0x12},
    };
    EXPECT_EQ(ReadI32LE(bytes.data()), 0x12345678);
}

TEST(CramInternal, TagPayloadHelpersSelectExpectedIntegerWidths)
{
    const std::vector<std::pair<std::int64_t, char>> cases{
        {-128, 'c'},  {255, 'C'},           {-32768, 's'},
        {65535, 'S'}, {-2147483648LL, 'i'}, {4000000000LL, 'I'},
    };

    for (const auto& [value, expectedType] : cases) {
        const auto encoded = EncodeTagValueToBamPayload(TagValue{value});
        EXPECT_EQ(encoded.Type, expectedType) << "value=" << value;
    }
}

TEST(CramInternal, TagPayloadHelpersRoundTripAllSupportedTypes)
{
    TagArray arrayValue{'C'};
    arrayValue.AppendUInt8(1);
    arrayValue.AppendUInt8(2);
    arrayValue.AppendUInt8(255);

    const std::vector<TagValue> values{
        TagValue{'Q'},
        TagValue{std::int64_t{-128}},
        TagValue{std::int64_t{255}},
        TagValue{std::int64_t{-32768}},
        TagValue{std::int64_t{65535}},
        TagValue{std::int64_t{-2147483648LL}},
        TagValue{std::int64_t{4000000000LL}},
        TagValue{3.25F},
        TagValue{std::string{"tag-text"}},
        TagValue{HexString{"0A0B"}},
        TagValue{arrayValue},
    };

    for (const auto& value : values) {
        const auto encoded = EncodeTagValueToBamPayload(value);
        const auto decoded = DecodeTagValueFromBamPayload(encoded.Type, encoded.Payload);
        if (const auto* floatValue = std::get_if<float>(&value)) {
            ASSERT_TRUE(std::holds_alternative<float>(decoded));
            EXPECT_FLOAT_EQ(std::get<float>(decoded), *floatValue);
            continue;
        }
        EXPECT_EQ(decoded, value);
    }
}

TEST(CramInternal, TagPayloadHelpersRejectInvalidPayloadWidths)
{
    const std::array<std::byte, 1> shortInt16{std::byte{0x01}};
    const std::array<std::byte, 3> shortInt32{std::byte{0x01}, std::byte{0x02}, std::byte{0x03}};
    const std::array<std::byte, 3> shortFloat{std::byte{0x00}, std::byte{0x00}, std::byte{0x80}};
    const std::array<std::byte, 6> malformedArray{std::byte{'C'},  std::byte{0x02},
                                                  std::byte{0x00}, std::byte{0x00},
                                                  std::byte{0x00}, std::byte{0x01}};

    EXPECT_THROW(DecodeTagValueFromBamPayload('s', shortInt16), std::runtime_error);
    EXPECT_THROW(DecodeTagValueFromBamPayload('i', shortInt32), std::runtime_error);
    EXPECT_THROW(DecodeTagValueFromBamPayload('f', shortFloat), std::runtime_error);
    EXPECT_THROW(DecodeTagValueFromBamPayload('B', malformedArray), std::runtime_error);
}

TEST(CramItf8, RoundTripSmallValues)
{
    for (const std::int32_t val : {0, 1, 42, 127}) {
        std::vector<std::byte> buf;
        WriteItf8(buf, val);
        EXPECT_EQ(std::size(buf), 1u) << "value=" << val;

        std::size_t bytesRead = 0;
        const auto decoded = ReadItf8(buf, bytesRead);
        EXPECT_EQ(decoded, val);
        EXPECT_EQ(bytesRead, 1u);
    }
}

TEST(CramItf8, RoundTripMediumValues)
{
    for (const std::int32_t val : {128, 255, 1000, 16383}) {
        std::vector<std::byte> buf;
        WriteItf8(buf, val);
        EXPECT_EQ(std::size(buf), 2u) << "value=" << val;

        std::size_t bytesRead = 0;
        const auto decoded = ReadItf8(buf, bytesRead);
        EXPECT_EQ(decoded, val);
        EXPECT_EQ(bytesRead, 2u);
    }
}

TEST(CramItf8, RoundTripLargeValues)
{
    for (const std::int32_t val : {16384, 100000, 2097151}) {
        std::vector<std::byte> buf;
        WriteItf8(buf, val);
        EXPECT_EQ(std::size(buf), 3u) << "value=" << val;

        std::size_t bytesRead = 0;
        const auto decoded = ReadItf8(buf, bytesRead);
        EXPECT_EQ(decoded, val);
        EXPECT_EQ(bytesRead, 3u);
    }
}

TEST(CramItf8, RoundTripNegativeValues)
{
    // Negative values are valid in ITF-8 (encoded as unsigned)
    for (const std::int32_t val : {-1, -100, -10000}) {
        std::vector<std::byte> buf;
        WriteItf8(buf, val);

        std::size_t bytesRead = 0;
        const auto decoded = ReadItf8(buf, bytesRead);
        EXPECT_EQ(decoded, val);
    }
}

TEST(CramLtf8, RoundTripSmallValues)
{
    for (const std::int64_t val : {0LL, 1LL, 42LL, 127LL}) {
        std::vector<std::byte> buf;
        WriteLtf8(buf, val);
        EXPECT_EQ(std::size(buf), 1u);

        std::size_t bytesRead = 0;
        const auto decoded = ReadLtf8(buf, bytesRead);
        EXPECT_EQ(decoded, val);
    }
}

TEST(CramLtf8, RoundTripLargeValues)
{
    for (const std::int64_t val : {128LL, 16384LL, 2097152LL, 268435456LL}) {
        std::vector<std::byte> buf;
        WriteLtf8(buf, val);

        std::size_t bytesRead = 0;
        const auto decoded = ReadLtf8(buf, bytesRead);
        EXPECT_EQ(decoded, val);
    }
}

// ===========================================================================
// File Definition Tests
// ===========================================================================

TEST(CramFileDefinition, RoundTrip)
{
    CramFileDefinition def;
    def.MajorVersion = 3;
    def.MinorVersion = 1;
    def.FileId[0] = std::byte{0x42};

    const auto bytes = SerializeFileDefinition(def);
    EXPECT_EQ(std::size(bytes), 26u);
    EXPECT_EQ(bytes[0], std::byte{0x43});  // 'C'
    EXPECT_EQ(bytes[1], std::byte{0x52});  // 'R'
    EXPECT_EQ(bytes[2], std::byte{0x41});  // 'A'
    EXPECT_EQ(bytes[3], std::byte{0x4d});  // 'M'
    EXPECT_EQ(bytes[4], std::byte{3});
    EXPECT_EQ(bytes[5], std::byte{1});

    const auto parsed = ParseFileDefinition(bytes);
    EXPECT_EQ(parsed.MajorVersion, 3u);
    EXPECT_EQ(parsed.MinorVersion, 1u);
    EXPECT_EQ(parsed.FileId[0], std::byte{0x42});
}

TEST(CramFileDefinition, InvalidMagicThrows)
{
    std::vector<std::byte> bad(26, std::byte{0});
    EXPECT_THROW(ParseFileDefinition(bad), std::runtime_error);
}

TEST(CramFileDefinition, TooShortThrows)
{
    std::vector<std::byte> small(10, std::byte{0});
    EXPECT_THROW(ParseFileDefinition(small), std::runtime_error);
}

// ===========================================================================
// EOF Marker Tests
// ===========================================================================

TEST(CramEof, RecognizeMarker)
{
    std::vector<std::byte> eof(std::begin(CRAM_EOF_MARKER), std::end(CRAM_EOF_MARKER));
    EXPECT_TRUE(IsCramEofMarker(eof));
}

TEST(CramEof, RejectTooShort)
{
    std::vector<std::byte> small(10, std::byte{0});
    EXPECT_FALSE(IsCramEofMarker(small));
}

TEST(CramEof, RejectBadData)
{
    std::vector<std::byte> bad(38, std::byte{0});
    EXPECT_FALSE(IsCramEofMarker(bad));
}

// ===========================================================================
// Block Tests
// ===========================================================================

TEST(CramBlock, SerializeAndParseRaw)
{
    CramBlock block;
    block.Method = CramBlockMethod::RAW;
    block.ContentType = CramBlockContentType::EXTERNAL_DATA;
    block.ContentId = 42;
    block.Data = {std::byte{0xDE}, std::byte{0xAD}, std::byte{0xBE}, std::byte{0xEF}};
    block.RawSize = 4;
    block.CompressedSize = 4;

    const auto serialized = SerializeBlock(block);

    std::size_t bytesRead = 0;
    const auto parsed = ParseBlock(serialized, bytesRead);
    EXPECT_EQ(bytesRead, std::size(serialized));
    EXPECT_EQ(parsed.Method, CramBlockMethod::RAW);
    EXPECT_EQ(parsed.ContentType, CramBlockContentType::EXTERNAL_DATA);
    EXPECT_EQ(parsed.ContentId, 42);
    EXPECT_EQ(parsed.CompressedSize, 4);
    EXPECT_EQ(parsed.RawSize, 4);
    EXPECT_EQ(std::size(parsed.Data), 4u);
    EXPECT_EQ(parsed.Data[0], std::byte{0xDE});
}

TEST(CramBlock, CrcMismatchThrows)
{
    CramBlock block;
    block.Method = CramBlockMethod::RAW;
    block.ContentType = CramBlockContentType::EXTERNAL_DATA;
    block.ContentId = 7;
    block.Data = {std::byte{0x10}, std::byte{0x20}, std::byte{0x30}};
    block.RawSize = 3;
    block.CompressedSize = 3;

    auto serialized = SerializeBlock(block);
    ASSERT_GT(std::size(serialized), 5u);
    serialized[4] ^= std::byte{0x01};  // corrupt payload/header byte without updating CRC

    std::size_t bytesRead = 0;
    EXPECT_THROW(ParseBlock(serialized, bytesRead), std::runtime_error);
}

TEST(CramBlock, NegativeCompressedSizeThrows)
{
    std::vector<std::byte> serialized;
    serialized.push_back(static_cast<std::byte>(CramBlockMethod::RAW));
    serialized.push_back(static_cast<std::byte>(CramBlockContentType::EXTERNAL_DATA));
    WriteItf8(serialized, 0);   // content id
    WriteItf8(serialized, -1);  // compressed size
    WriteItf8(serialized, 0);   // raw size

    std::size_t bytesRead = 0;
    EXPECT_THROW(ParseBlock(serialized, bytesRead), std::runtime_error);
}

// ===========================================================================
// Container Header Tests
// ===========================================================================

TEST(CramContainerHeader, SerializeAndParse)
{
    CramContainerHeader hdr;
    hdr.Length = 100;
    hdr.RefSeqId = 0;
    hdr.StartPos = 1;
    hdr.AlignmentSpan = 1000;
    hdr.NumRecords = 50;
    hdr.RecordCounter = 0;
    hdr.Bases = 5000;
    hdr.NumBlocks = 3;
    hdr.Landmarks = {10, 50};

    const auto serialized = SerializeContainerHeader(hdr);

    std::size_t bytesRead = 0;
    const auto parsed = ParseContainerHeader(serialized, bytesRead);
    EXPECT_EQ(bytesRead, std::size(serialized));
    EXPECT_EQ(parsed.Length, 100);
    EXPECT_EQ(parsed.RefSeqId, 0);
    EXPECT_EQ(parsed.StartPos, 1);
    EXPECT_EQ(parsed.AlignmentSpan, 1000);
    EXPECT_EQ(parsed.NumRecords, 50);
    EXPECT_EQ(parsed.RecordCounter, 0);
    EXPECT_EQ(parsed.Bases, 5000);
    EXPECT_EQ(parsed.NumBlocks, 3);
    ASSERT_EQ(std::size(parsed.Landmarks), 2u);
    EXPECT_EQ(parsed.Landmarks[0], 10);
    EXPECT_EQ(parsed.Landmarks[1], 50);
}

TEST(CramContainerHeader, CrcMismatchThrows)
{
    CramContainerHeader hdr;
    hdr.Length = 11;
    hdr.RefSeqId = 0;
    hdr.StartPos = 1;
    hdr.AlignmentSpan = 10;
    hdr.NumRecords = 1;
    hdr.RecordCounter = 0;
    hdr.Bases = 4;
    hdr.NumBlocks = 1;
    hdr.Landmarks = {0};

    auto serialized = SerializeContainerHeader(hdr);
    ASSERT_GT(std::size(serialized), 4u);
    serialized[4] ^= std::byte{0x01};  // corrupt payload without updating CRC

    std::size_t bytesRead = 0;
    EXPECT_THROW(ParseContainerHeader(serialized, bytesRead), std::runtime_error);
}

TEST(CramContainerHeader, NegativeLandmarkCountThrows)
{
    std::vector<std::byte> serialized(4, std::byte{0});
    std::int32_t length{0};
    std::memcpy(serialized.data(), &length, 4);

    WriteItf8(serialized, 0);   // ref
    WriteItf8(serialized, 1);   // start
    WriteItf8(serialized, 10);  // span
    WriteItf8(serialized, 0);   // num_records
    WriteLtf8(serialized, 0);   // record_counter
    WriteLtf8(serialized, 0);   // bases
    WriteItf8(serialized, 1);   // num_blocks
    WriteItf8(serialized, -1);  // landmark count (invalid)

    const auto crc = libdeflate_crc32(0, serialized.data(), std::size(serialized));
    const auto crcPos = std::size(serialized);
    serialized.resize(crcPos + 4);
    std::memcpy(serialized.data() + crcPos, &crc, 4);

    std::size_t bytesRead = 0;
    EXPECT_THROW(ParseContainerHeader(serialized, bytesRead), std::runtime_error);
}

// ===========================================================================
// Slice Header Tests
// ===========================================================================

TEST(CramSliceHeader, SerializeAndParse)
{
    CramSliceHeader hdr;
    hdr.RefSeqId = 0;
    hdr.AlignmentStart = 100;
    hdr.AlignmentSpan = 500;
    hdr.NumRecords = 25;
    hdr.RecordCounter = 0;
    hdr.NumBlocks = 3;
    hdr.BlockContentIds = {1, 2, 3};
    hdr.EmbeddedRefBlockId = -1;

    const auto serialized = SerializeSliceHeader(hdr);
    const auto parsed = ParseSliceHeader(serialized);

    EXPECT_EQ(parsed.RefSeqId, 0);
    EXPECT_EQ(parsed.AlignmentStart, 100);
    EXPECT_EQ(parsed.AlignmentSpan, 500);
    EXPECT_EQ(parsed.NumRecords, 25);
    EXPECT_EQ(parsed.NumBlocks, 3);
    ASSERT_EQ(std::size(parsed.BlockContentIds), 3u);
    EXPECT_EQ(parsed.EmbeddedRefBlockId, -1);
}

TEST(CramSliceHeader, OptionalTagsRoundTrip)
{
    CramSliceHeader hdr;
    hdr.RefSeqId = 1;
    hdr.AlignmentStart = 50;
    hdr.AlignmentSpan = 200;
    hdr.NumRecords = 2;
    hdr.RecordCounter = 10;
    hdr.NumBlocks = 2;
    hdr.BlockContentIds = {100, 101};
    hdr.EmbeddedRefBlockId = -1;
    hdr.OptionalTags = {std::byte{'N'}, std::byte{'M'}, std::byte{'i'}, std::byte{5},
                        std::byte{0},   std::byte{0},   std::byte{0}};

    const auto serialized = SerializeSliceHeader(hdr);
    const auto parsed = ParseSliceHeader(serialized);

    EXPECT_EQ(parsed.RefSeqId, hdr.RefSeqId);
    EXPECT_EQ(parsed.AlignmentStart, hdr.AlignmentStart);
    EXPECT_EQ(parsed.AlignmentSpan, hdr.AlignmentSpan);
    EXPECT_EQ(parsed.NumRecords, hdr.NumRecords);
    EXPECT_EQ(parsed.RecordCounter, hdr.RecordCounter);
    EXPECT_EQ(parsed.NumBlocks, hdr.NumBlocks);
    EXPECT_EQ(parsed.BlockContentIds, hdr.BlockContentIds);
    EXPECT_EQ(parsed.EmbeddedRefBlockId, hdr.EmbeddedRefBlockId);
    EXPECT_EQ(parsed.OptionalTags, hdr.OptionalTags);
}

TEST(CramSlice, DefaultConstructionIsEmpty)
{
    CramSlice slice;
    EXPECT_EQ(slice.Header.NumRecords, 0);
    EXPECT_EQ(slice.CoreBlock.ContentType, CramBlockContentType::FILE_HEADER);
    EXPECT_TRUE(std::empty(slice.ExternalBlocks));
}

TEST(CramContainer, DefaultConstructionIsEmpty)
{
    CramContainer container;
    EXPECT_EQ(container.Header.NumRecords, 0);
    EXPECT_TRUE(std::empty(container.CompressionHeader.DataSeriesEncodings));
    EXPECT_TRUE(std::empty(container.Slices));
}

TEST(CramSlice, SerializeRoundTrip)
{
    CramSlice slice;
    slice.Header.RefSeqId = 0;
    slice.Header.AlignmentStart = 100;
    slice.Header.AlignmentSpan = 50;
    slice.Header.NumRecords = 2;
    slice.Header.RecordCounter = 0;
    slice.Header.NumBlocks = 2;
    slice.Header.BlockContentIds = {1};

    slice.CoreBlock.Method = CramBlockMethod::RAW;
    slice.CoreBlock.ContentType = CramBlockContentType::CORE_DATA;
    slice.CoreBlock.ContentId = 0;
    slice.CoreBlock.RawSize = 0;
    slice.CoreBlock.CompressedSize = 0;

    CramBlock extBlock;
    extBlock.Method = CramBlockMethod::RAW;
    extBlock.ContentType = CramBlockContentType::EXTERNAL_DATA;
    extBlock.ContentId = 1;
    extBlock.Data = {std::byte{0xAA}, std::byte{0xBB}};
    extBlock.RawSize = 2;
    extBlock.CompressedSize = 2;
    slice.ExternalBlocks.push_back(std::move(extBlock));

    const auto bytes = SerializeSlice(slice);
    ASSERT_FALSE(std::empty(bytes));

    std::size_t bytesRead = 0;
    const auto headerBlock = ParseBlock(bytes, bytesRead);
    EXPECT_EQ(headerBlock.ContentType, CramBlockContentType::SLICE_HEADER);

    const auto parsedHeader = ParseSliceHeader(headerBlock.Data);
    EXPECT_EQ(parsedHeader.RefSeqId, 0);
    EXPECT_EQ(parsedHeader.AlignmentStart, 100);
    EXPECT_EQ(parsedHeader.NumRecords, 2);
    EXPECT_EQ(SerializedSliceSize(slice), std::size(bytes));
}

TEST(CramContainer, SerializeAndParseRoundTrip)
{
    CramCompressionHeader compressionHeader;
    compressionHeader.PreservationMap.ReadNamesIncluded = true;
    compressionHeader.PreservationMap.ApDelta = false;
    compressionHeader.PreservationMap.ReferenceRequired = false;

    CramEncodingDescriptor extDesc;
    extDesc.CodecId = CramCodecId::EXTERNAL;
    WriteItf8(extDesc.Parameters, 1);
    compressionHeader.DataSeriesEncodings.emplace_back(CramDataSeries::BF, extDesc);

    CramContainer container;
    container.Header.RefSeqId = -1;
    container.Header.StartPos = 1;
    container.Header.AlignmentSpan = 0;
    container.Header.NumRecords = 2;
    container.Header.RecordCounter = 0;
    container.Header.Bases = 0;
    container.CompressionHeader = std::move(compressionHeader);
    container.Slices = {MakeRawExternalSlice(0, 1), MakeRawExternalSlice(1, 2)};

    const auto serialized = SerializeContainer(container);
    ASSERT_FALSE(std::empty(serialized));

    std::size_t headerBytesRead = 0;
    const auto parsedHeader = ParseContainerHeader(serialized, headerBytesRead);
    const auto parsed = ParseContainer(
        parsedHeader, std::span<const std::byte>{serialized}.subspan(headerBytesRead));

    ASSERT_EQ(std::size(parsed.Slices), 2u);
    EXPECT_EQ(parsed.Slices[0].Header.RecordCounter, 0);
    EXPECT_EQ(parsed.Slices[1].Header.RecordCounter, 1);
    EXPECT_EQ(parsed.Slices[0].CoreBlock.ContentType, CramBlockContentType::CORE_DATA);
    EXPECT_EQ(parsed.Slices[1].ExternalBlocks.front().ContentId, 2);
}

// ===========================================================================
// Compression Header Tests
// ===========================================================================

TEST(CramCompressionHeader, SerializeAndParse)
{
    CramCompressionHeader hdr;
    hdr.PreservationMap.ReadNamesIncluded = true;
    hdr.PreservationMap.ApDelta = false;
    hdr.PreservationMap.ReferenceRequired = false;

    // Add one data series encoding
    CramEncodingDescriptor extDesc;
    extDesc.CodecId = CramCodecId::EXTERNAL;
    WriteItf8(extDesc.Parameters, 1);
    hdr.DataSeriesEncodings.emplace_back(CramDataSeries::BF, std::move(extDesc));

    const auto serialized = SerializeCompressionHeader(hdr);
    const auto parsed = ParseCompressionHeader(serialized);

    EXPECT_EQ(parsed.PreservationMap.ReadNamesIncluded, true);
    EXPECT_EQ(parsed.PreservationMap.ApDelta, false);
    EXPECT_EQ(parsed.PreservationMap.ReferenceRequired, false);
    ASSERT_EQ(std::size(parsed.DataSeriesEncodings), 1u);
    EXPECT_EQ(parsed.DataSeriesEncodings[0].first, CramDataSeries::BF);
    EXPECT_EQ(parsed.DataSeriesEncodings[0].second.CodecId, CramCodecId::EXTERNAL);
}

TEST(CramCompressionHeader, TruncatedDataThrows)
{
    const CramCompressionHeader hdr;
    const auto serialized = SerializeCompressionHeader(hdr);
    ASSERT_FALSE(std::empty(serialized));

    auto truncated = serialized;
    truncated.pop_back();
    EXPECT_THROW(ParseCompressionHeader(truncated), std::runtime_error);
}

// ===========================================================================
// Bit Reader / Writer Tests
// ===========================================================================

TEST(CramBit, ReadWriteSingleBits)
{
    CramBitWriter writer;
    writer.WriteBit(1);
    writer.WriteBit(0);
    writer.WriteBit(1);
    writer.WriteBit(1);
    writer.Flush();

    const auto data = writer.Data();
    CramBitReader reader{data};
    EXPECT_EQ(reader.ReadBit(), 1);
    EXPECT_EQ(reader.ReadBit(), 0);
    EXPECT_EQ(reader.ReadBit(), 1);
    EXPECT_EQ(reader.ReadBit(), 1);
}

TEST(CramBit, ReadWriteMultiBits)
{
    CramBitWriter writer;
    writer.WriteBits(42, 8);
    writer.WriteBits(7, 4);
    writer.Flush();

    const auto data = writer.Data();
    CramBitReader reader{data};
    EXPECT_EQ(reader.ReadBits(8), 42);
    EXPECT_EQ(reader.ReadBits(4), 7);
}

// ===========================================================================
// External Block Store Tests
// ===========================================================================

TEST(CramExternalBlockStore, WriteAndRead)
{
    CramExternalBlockStore store;

    std::vector<std::byte> data = {std::byte{0x01}, std::byte{0x02}, std::byte{0x03}};
    store.AddBlock(1, data);

    auto result = store.ReadBytes(1, 2);
    ASSERT_EQ(std::size(result), 2u);
    EXPECT_EQ(result[0], std::byte{0x01});
    EXPECT_EQ(result[1], std::byte{0x02});

    auto result2 = store.ReadBytes(1, 1);
    EXPECT_EQ(result2[0], std::byte{0x03});
}

TEST(CramExternalBlockStore, WriteItf8AndRead)
{
    CramExternalBlockStore store;
    store.WriteItf8(1, 42);
    store.WriteItf8(1, 1000);

    store.ResetPositions();

    EXPECT_EQ(store.ReadItf8(1), 42);
    EXPECT_EQ(store.ReadItf8(1), 1000);
}

TEST(CramExternalBlockStore, ContentIds)
{
    CramExternalBlockStore store;
    store.WriteItf8(3, 0);
    store.WriteItf8(1, 0);
    store.WriteItf8(5, 0);

    auto ids = store.ContentIds();
    ASSERT_EQ(std::size(ids), 3u);
    EXPECT_EQ(ids[0], 1);
    EXPECT_EQ(ids[1], 3);
    EXPECT_EQ(ids[2], 5);
}

// ===========================================================================
// Codec Tests
// ===========================================================================

TEST(CramCodec, ExternalCodecRoundTrip)
{
    CramExternalBlockStore store;
    CramBitReader dummyReader{{}};
    CramBitWriter dummyWriter;

    ExternalCodec codec{1};

    codec.EncodeInt(42, dummyWriter, store);
    codec.EncodeInt(-1, dummyWriter, store);

    store.ResetPositions();

    EXPECT_EQ(codec.DecodeInt(dummyReader, store), 42);
    EXPECT_EQ(codec.DecodeInt(dummyReader, store), -1);
}

TEST(CramCodec, ExternalCodecByteRoundTrip)
{
    CramExternalBlockStore store;
    CramBitReader dummyReader{{}};
    CramBitWriter dummyWriter;

    ExternalCodec codec{1};
    codec.EncodeByte(std::byte{0xAB}, dummyWriter, store);

    store.ResetPositions();

    EXPECT_EQ(codec.DecodeByte(dummyReader, store), std::byte{0xAB});
}

TEST(CramCodec, ExternalCodecByteArrayRoundTrip)
{
    CramExternalBlockStore store;
    CramBitReader dummyReader{{}};
    CramBitWriter dummyWriter;

    ExternalCodec codec{1};
    std::vector<std::byte> original = {std::byte{'H'}, std::byte{'i'}};
    codec.EncodeByteArray(original, dummyWriter, store);

    store.ResetPositions();

    auto decoded = codec.DecodeByteArray(dummyReader, store);
    ASSERT_EQ(std::size(decoded), 2u);
    EXPECT_EQ(decoded[0], std::byte{'H'});
    EXPECT_EQ(decoded[1], std::byte{'i'});
}

TEST(CramCodec, DecodeKindMatchesCodecCapabilities)
{
    const NullCodec nullCodec;
    const ExternalCodec externalCodec{1};
    const HuffmanCodec huffmanCodec({42}, {0});
    const BetaCodec betaCodec{0, 8};
    const GammaCodec gammaCodec{0};
    const SubexpCodec subexpCodec{0, 2};
    const ByteArrayLenCodec byteArrayLenCodec{std::make_unique<ExternalCodec>(1),
                                              std::make_unique<ExternalCodec>(1)};
    const ByteArrayStopCodec byteArrayStopCodec{std::byte{0}, 1};

    EXPECT_EQ(nullCodec.DecodeKind(), CramCodecDecodeKind::SCALAR);
    EXPECT_EQ(externalCodec.DecodeKind(), CramCodecDecodeKind::SCALAR);
    EXPECT_EQ(huffmanCodec.DecodeKind(), CramCodecDecodeKind::SCALAR);
    EXPECT_EQ(betaCodec.DecodeKind(), CramCodecDecodeKind::SCALAR);
    EXPECT_EQ(gammaCodec.DecodeKind(), CramCodecDecodeKind::SCALAR);
    EXPECT_EQ(subexpCodec.DecodeKind(), CramCodecDecodeKind::SCALAR);
    EXPECT_EQ(byteArrayLenCodec.DecodeKind(), CramCodecDecodeKind::BYTE_ARRAY);
    EXPECT_EQ(byteArrayStopCodec.DecodeKind(), CramCodecDecodeKind::BYTE_ARRAY);
}

TEST(CramCodec, SubexpRoundTripCanonical)
{
    // Encode then decode a series; sequential decoding only stays aligned if each value
    // consumes exactly the bits written (the over-read bug misaligns subsequent values).
    // This is the canonical Golomb-subexponential bitstream htslib decodes.
    const std::vector<std::int32_t> values{0, 1, 2, 3, 4, 6, 7, 8, 15, 16, 31, 100, 255, 500};
    for (const std::int32_t k : {0, 1, 2, 3, 4}) {
        for (const std::int32_t offset : {0, 3}) {
            SubexpCodec codec{offset, k};
            CramBitWriter writer;
            CramExternalBlockStore store;
            for (const std::int32_t v : values) {
                codec.EncodeInt(v, writer, store);
            }
            writer.Flush();

            const auto data = writer.Data();
            CramBitReader reader{data};
            for (const std::int32_t v : values) {
                EXPECT_EQ(codec.DecodeInt(reader, store), v) << "k=" << k << " offset=" << offset;
            }
        }
    }
}

TEST(CramCodec, HuffmanSingleSymbol)
{
    HuffmanCodec codec({42}, {0});

    std::vector<std::byte> emptyData;
    CramBitReader reader{emptyData};
    CramExternalBlockStore store;

    // Single symbol → always returns that symbol
    EXPECT_EQ(codec.DecodeInt(reader, store), 42);
}

TEST(CramCodec, HuffmanMixedLengthRoundTrip)
{
    HuffmanCodec codec({10, 20, 30, 40}, {1, 3, 3, 2});

    CramBitWriter writer;
    CramExternalBlockStore store;
    for (const std::int32_t value : {40, 10, 30, 20, 10}) {
        codec.EncodeInt(value, writer, store);
    }
    writer.Flush();

    const auto data = writer.Data();
    CramBitReader reader{data};
    EXPECT_EQ(codec.DecodeInt(reader, store), 40);
    EXPECT_EQ(codec.DecodeInt(reader, store), 10);
    EXPECT_EQ(codec.DecodeInt(reader, store), 30);
    EXPECT_EQ(codec.DecodeInt(reader, store), 20);
    EXPECT_EQ(codec.DecodeInt(reader, store), 10);
}

TEST(CramCodec, HuffmanMalformedCodeThrowsDecodeFailure)
{
    HuffmanCodec codec({1, 2}, {1, 2});

    CramBitWriter writer;
    for (int i = 0; i < 33; ++i) {
        writer.WriteBit(1);
    }
    writer.Flush();

    const auto bits = writer.Data();
    CramBitReader reader{bits};
    CramExternalBlockStore store;

    try {
        static_cast<void>(codec.DecodeInt(reader, store));
        FAIL() << "expected malformed code decode failure";
    } catch (const std::runtime_error& err) {
        EXPECT_EQ(std::string{err.what()}, "HuffmanCodec: decode failed — no matching code");
    }
}

TEST(CramCodec, HuffmanDecodeUsesIndexedLookupWithoutLinearEntryScan)
{
    std::vector<std::byte> sourceBytes;
    for (const auto& path : {std::filesystem::path{"../src/CramCodec.cpp"},
                             std::filesystem::path{"src/CramCodec.cpp"}}) {
        if (std::filesystem::exists(path)) {
            sourceBytes = ReadFileBytes(path);
            break;
        }
    }
    ASSERT_FALSE(sourceBytes.empty()) << "failed to read CramCodec.cpp source";

    const std::string source(reinterpret_cast<const char*>(sourceBytes.data()), sourceBytes.size());
    const std::size_t decodeStart = source.find("std::int32_t HuffmanCodec::DecodeInt");
    const std::size_t decodeEnd = source.find("std::byte HuffmanCodec::DecodeByte", decodeStart);
    ASSERT_NE(decodeStart, std::string::npos);
    ASSERT_NE(decodeEnd, std::string::npos);

    const std::string decodeBody = source.substr(decodeStart, decodeEnd - decodeStart);
    EXPECT_EQ(decodeBody.find("for (const auto& e : entries_)"), std::string::npos);
}

TEST(CramCodec, DecodeTagPayloadUsesExplicitCodecKindWithoutExceptionProbing)
{
    std::vector<std::byte> sourceBytes;
    for (const auto& path : {std::filesystem::path{"../src/CramReader.cpp"},
                             std::filesystem::path{"src/CramReader.cpp"}}) {
        if (std::filesystem::exists(path)) {
            sourceBytes = ReadFileBytes(path);
            break;
        }
    }
    ASSERT_FALSE(sourceBytes.empty()) << "failed to read CramReader.cpp source";

    const std::string source(reinterpret_cast<const char*>(sourceBytes.data()), sourceBytes.size());
    const std::size_t decodeStart = source.find("std::vector<std::byte> DecodeTagPayload");
    const std::size_t decodeEnd = source.find("std::array<std::array<char, 4>, 5>", decodeStart);
    ASSERT_NE(decodeStart, std::string::npos);
    ASSERT_NE(decodeEnd, std::string::npos);

    const std::string decodeBody = source.substr(decodeStart, decodeEnd - decodeStart);
    EXPECT_NE(decodeBody.find("DecodeKind()"), std::string::npos);
    EXPECT_EQ(decodeBody.find("tryDecodeByteArray"), std::string::npos);
    EXPECT_EQ(decodeBody.find("catch (...)"), std::string::npos);
}

TEST(CramCodec, BetaCodecRoundTrip)
{
    BetaCodec codec{0, 8};  // 8-bit values, offset 0

    CramBitWriter writer;
    CramExternalBlockStore store;

    codec.EncodeInt(42, writer, store);
    codec.EncodeInt(255, writer, store);
    writer.Flush();

    const auto data = writer.Data();
    CramBitReader reader{data};
    EXPECT_EQ(codec.DecodeInt(reader, store), 42);
    EXPECT_EQ(codec.DecodeInt(reader, store), 255);
}

TEST(CramCodec, GammaCodecRoundTrip)
{
    GammaCodec codec{0};

    CramBitWriter writer;
    CramExternalBlockStore store;

    codec.EncodeInt(1, writer, store);
    codec.EncodeInt(5, writer, store);
    codec.EncodeInt(100, writer, store);
    writer.Flush();

    const auto data = writer.Data();
    CramBitReader reader{data};
    EXPECT_EQ(codec.DecodeInt(reader, store), 1);
    EXPECT_EQ(codec.DecodeInt(reader, store), 5);
    EXPECT_EQ(codec.DecodeInt(reader, store), 100);
}

TEST(CramCodec, NullCodecReturnsDefaults)
{
    NullCodec codec;

    std::vector<std::byte> emptyData;
    CramBitReader reader{emptyData};
    CramExternalBlockStore store;

    EXPECT_EQ(codec.DecodeInt(reader, store), 0);
    EXPECT_EQ(codec.DecodeByte(reader, store), std::byte{0});
    EXPECT_TRUE(std::empty(codec.DecodeByteArray(reader, store)));
}

TEST(CramCodec, ByteArrayStopCodecRoundTrip)
{
    CramExternalBlockStore store;
    CramBitReader dummyReader{{}};
    CramBitWriter dummyWriter;

    ByteArrayStopCodec codec{std::byte{0x00}, 1};

    std::vector<std::byte> original = {std::byte{'A'}, std::byte{'B'}, std::byte{'C'}};
    codec.EncodeByteArray(original, dummyWriter, store);

    store.ResetPositions();

    auto decoded = codec.DecodeByteArray(dummyReader, store);
    ASSERT_EQ(std::size(decoded), 3u);
    EXPECT_EQ(decoded[0], std::byte{'A'});
    EXPECT_EQ(decoded[1], std::byte{'B'});
    EXPECT_EQ(decoded[2], std::byte{'C'});
}

TEST(CramCodec, CreateCodecFromDescriptor)
{
    // External codec
    CramEncodingDescriptor desc;
    desc.CodecId = CramCodecId::EXTERNAL;
    WriteItf8(desc.Parameters, 5);

    auto codec = CreateCodec(desc);
    ASSERT_TRUE(codec);

    // Null codec
    CramEncodingDescriptor nullDesc;
    nullDesc.CodecId = CramCodecId::NONE;
    auto nullCodec = CreateCodec(nullDesc);
    ASSERT_TRUE(nullCodec);
}

// ===========================================================================
// Compression Tests
// ===========================================================================

TEST(CramCompression, GzipRoundTrip)
{
    // Use a longer, repetitive string to ensure gzip actually compresses
    std::string original;
    for (int i = 0; i < 100; ++i) {
        original += "Hello, CRAM compression world! This is a test string for gzip. ";
    }
    const std::span<const std::byte> input{reinterpret_cast<const std::byte*>(original.data()),
                                           std::size(original)};

    const auto compressed = CramGzipCompress(input);
    EXPECT_LT(std::size(compressed), std::size(original));

    const auto decompressed = CramGzipDecompress(compressed, std::size(original));
    ASSERT_EQ(std::size(decompressed), std::size(original));

    const std::string result(reinterpret_cast<const char*>(decompressed.data()),
                             std::size(decompressed));
    EXPECT_EQ(result, original);
}

TEST(CramCompression, RawPassthrough)
{
    auto& registry = CompressionRegistry::Instance();
    std::vector<std::byte> data = {std::byte{1}, std::byte{2}, std::byte{3}};

    auto compressed = registry.Compress(CramBlockMethod::RAW, data);
    EXPECT_EQ(compressed, data);

    auto decompressed = registry.Decompress(CramBlockMethod::RAW, data, 3);
    EXPECT_EQ(decompressed, data);
}

TEST(CramCompression, RegistryHasBuiltins)
{
    auto& registry = CompressionRegistry::Instance();
    EXPECT_TRUE(registry.HasMethod(0));  // raw
    EXPECT_TRUE(registry.HasMethod(1));  // gzip
    EXPECT_TRUE(registry.HasMethod(4));  // rANS4x8
    EXPECT_TRUE(registry.HasMethod(5));  // rANS4x16
    EXPECT_TRUE(registry.HasMethod(6));  // adaptive arithmetic
    EXPECT_TRUE(registry.HasMethod(7));  // fqzcomp
    EXPECT_TRUE(registry.HasMethod(8));  // name tokeniser
#ifdef PBSAMOA_HAVE_BZIP2
    EXPECT_TRUE(registry.HasMethod(2));  // bzip2
#endif
#ifdef PBSAMOA_HAVE_LZMA
    EXPECT_TRUE(registry.HasMethod(3));  // lzma
#endif
}

TEST(CramCompression, CannotOverrideBuiltins)
{
    auto& registry = CompressionRegistry::Instance();
    EXPECT_THROW(registry.Register(0, EmptyCompressionResult, EmptyDecompressionResult),
                 std::runtime_error);
    EXPECT_THROW(registry.Register(8, EmptyCompressionResult, EmptyDecompressionResult),
                 std::runtime_error);
}

TEST(CramCompression, CustomRegistration)
{
    auto& registry = CompressionRegistry::Instance();

    // Register a custom "XOR" compression as method 200
    registry.Register(200, XorCompress, XorDecompress);

    EXPECT_TRUE(registry.HasMethod(200));

    std::vector<std::byte> original = {std::byte{0x12}, std::byte{0x34}};
    auto compressed = registry.Compress(static_cast<CramBlockMethod>(200), original);
    auto decompressed =
        registry.Decompress(static_cast<CramBlockMethod>(200), compressed, std::size(original));
    EXPECT_EQ(decompressed, original);
}

TEST(CramCompression, BlockCompressDecompress)
{
    CramBlock block;
    block.ContentType = CramBlockContentType::EXTERNAL_DATA;
    block.Data = {std::byte{0x01}, std::byte{0x02}, std::byte{0x03}, std::byte{0x04},
                  std::byte{0x05}, std::byte{0x06}, std::byte{0x07}, std::byte{0x08}};
    block.RawSize = static_cast<std::int32_t>(std::size(block.Data));

    CompressCramBlock(block, CramBlockMethod::GZIP);
    EXPECT_EQ(block.Method, CramBlockMethod::GZIP);

    DecompressCramBlock(block);
    EXPECT_EQ(block.Method, CramBlockMethod::RAW);
    EXPECT_EQ(std::size(block.Data), 8u);
    EXPECT_EQ(block.Data[0], std::byte{0x01});
    EXPECT_EQ(block.Data[7], std::byte{0x08});
}

TEST(CramCompression, BlockCompressDecompressWithReusableGzipContexts)
{
    const LibdeflateCompressorPtr compressor{libdeflate_alloc_compressor(6)};
    ASSERT_TRUE(compressor);
    const LibdeflateDecompressorPtr decompressor{libdeflate_alloc_decompressor()};
    ASSERT_TRUE(decompressor);

    CramBlock block;
    block.ContentType = CramBlockContentType::EXTERNAL_DATA;
    block.Data = {std::byte{0x11}, std::byte{0x22}, std::byte{0x33}, std::byte{0x44},
                  std::byte{0x55}, std::byte{0x66}, std::byte{0x77}, std::byte{0x88}};
    block.RawSize = static_cast<std::int32_t>(std::size(block.Data));

    CompressCramBlock(block, CramBlockMethod::GZIP, compressor.get());
    EXPECT_EQ(block.Method, CramBlockMethod::GZIP);

    DecompressCramBlock(block, decompressor.get());
    EXPECT_EQ(block.Method, CramBlockMethod::RAW);
    EXPECT_EQ(std::size(block.Data), 8u);
    EXPECT_EQ(block.Data[0], std::byte{0x11});
    EXPECT_EQ(block.Data[7], std::byte{0x88});
}

TEST(CramReaderInfrastructure, DecodeRecordUsesContextStruct)
{
    std::vector<std::byte> sourceBytes;
    for (const auto& path : {std::filesystem::path{"../src/CramReader.cpp"},
                             std::filesystem::path{"src/CramReader.cpp"}}) {
        if (std::filesystem::exists(path)) {
            sourceBytes = ReadFileBytes(path);
            break;
        }
    }
    ASSERT_FALSE(sourceBytes.empty()) << "failed to read CramReader.cpp source";

    const std::string source(reinterpret_cast<const char*>(sourceBytes.data()), sourceBytes.size());
    EXPECT_NE(source.find("struct DecodeRecordContext"), std::string::npos);
    EXPECT_NE(source.find("DecodeRecord(BamRecord& record"), std::string::npos);
    EXPECT_NE(source.find("DecodeRecordContext& ctx"), std::string::npos);
}

TEST(CramCompression, BlockCompressDecompressRans4x8)
{
    CramBlock block;
    block.ContentType = CramBlockContentType::EXTERNAL_DATA;
    block.Data = {std::byte{'A'}, std::byte{'C'}, std::byte{'G'}, std::byte{'T'},
                  std::byte{'A'}, std::byte{'C'}, std::byte{'G'}, std::byte{'T'}};
    block.RawSize = static_cast<std::int32_t>(std::size(block.Data));

    CompressCramBlock(block, CramBlockMethod::RANS4X8);
    EXPECT_EQ(block.Method, CramBlockMethod::RANS4X8);

    DecompressCramBlock(block);
    EXPECT_EQ(block.Method, CramBlockMethod::RAW);
    EXPECT_EQ(std::size(block.Data), 8u);
}

TEST(CramCompression, BlockCompressDecompressRans4x16)
{
    CramBlock block;
    block.ContentType = CramBlockContentType::EXTERNAL_DATA;
    block.Data = {std::byte{'A'}, std::byte{'A'}, std::byte{'A'}, std::byte{'C'},
                  std::byte{'C'}, std::byte{'C'}, std::byte{'G'}, std::byte{'G'}};
    block.RawSize = static_cast<std::int32_t>(std::size(block.Data));

    CompressCramBlock(block, CramBlockMethod::RANS4X16);
    EXPECT_EQ(block.Method, CramBlockMethod::RANS4X16);

    DecompressCramBlock(block);
    EXPECT_EQ(block.Method, CramBlockMethod::RAW);
    EXPECT_EQ(std::size(block.Data), 8u);
}

TEST(CramCompression, BlockCompressDecompressAdaptiveArith)
{
    CramBlock block;
    block.ContentType = CramBlockContentType::EXTERNAL_DATA;
    block.Data = {std::byte{'Q'}, std::byte{'Q'}, std::byte{'Q'}, std::byte{'Q'},
                  std::byte{'R'}, std::byte{'R'}, std::byte{'S'}, std::byte{'S'}};
    block.RawSize = static_cast<std::int32_t>(std::size(block.Data));

    CompressCramBlock(block, CramBlockMethod::ADAPTIVE_ARITH);
    EXPECT_EQ(block.Method, CramBlockMethod::ADAPTIVE_ARITH);

    DecompressCramBlock(block);
    EXPECT_EQ(block.Method, CramBlockMethod::RAW);
    EXPECT_EQ(std::size(block.Data), 8u);
}

TEST(CramCompression, BlockCompressDecompressFqzcomp)
{
    CramBlock block;
    block.ContentType = CramBlockContentType::EXTERNAL_DATA;
    block.Data = {std::byte{30}, std::byte{31}, std::byte{32}, std::byte{33},
                  std::byte{34}, std::byte{35}, std::byte{36}, std::byte{37}};
    block.RawSize = static_cast<std::int32_t>(std::size(block.Data));

    CompressCramBlock(block, CramBlockMethod::FQZCOMP);
    EXPECT_EQ(block.Method, CramBlockMethod::FQZCOMP);

    DecompressCramBlock(block);
    EXPECT_EQ(block.Method, CramBlockMethod::RAW);
    EXPECT_EQ(std::size(block.Data), 8u);
}

TEST(CramCompression, BlockCompressDecompressNameTokeniser)
{
    CramBlock block;
    block.ContentType = CramBlockContentType::EXTERNAL_DATA;
    block.Data = {std::byte{'r'}, std::byte{'e'}, std::byte{'a'}, std::byte{'d'},
                  std::byte{'1'}, std::byte{0},   std::byte{'r'}, std::byte{'e'},
                  std::byte{'a'}, std::byte{'d'}, std::byte{'2'}, std::byte{0}};
    block.RawSize = static_cast<std::int32_t>(std::size(block.Data));

    CompressCramBlock(block, CramBlockMethod::NAME_TOKENISER);
    EXPECT_EQ(block.Method, CramBlockMethod::NAME_TOKENISER);

    DecompressCramBlock(block);
    EXPECT_EQ(block.Method, CramBlockMethod::RAW);
    EXPECT_EQ(std::size(block.Data), 12u);
}

#ifdef PBSAMOA_HAVE_BZIP2
TEST(CramCompression, BlockCompressDecompressBzip2)
{
    CramBlock block;
    block.ContentType = CramBlockContentType::EXTERNAL_DATA;
    block.Data = {std::byte{0x10}, std::byte{0x11}, std::byte{0x12}, std::byte{0x13},
                  std::byte{0x14}, std::byte{0x15}, std::byte{0x16}, std::byte{0x17}};
    block.RawSize = static_cast<std::int32_t>(std::size(block.Data));

    CompressCramBlock(block, CramBlockMethod::BZIP2);
    EXPECT_EQ(block.Method, CramBlockMethod::BZIP2);

    DecompressCramBlock(block);
    EXPECT_EQ(block.Method, CramBlockMethod::RAW);
    EXPECT_EQ(std::size(block.Data), 8u);
    EXPECT_EQ(block.Data[0], std::byte{0x10});
    EXPECT_EQ(block.Data[7], std::byte{0x17});
}
#endif

#ifdef PBSAMOA_HAVE_LZMA
TEST(CramCompression, BlockCompressDecompressLzma)
{
    CramBlock block;
    block.ContentType = CramBlockContentType::EXTERNAL_DATA;
    block.Data = {std::byte{0x21}, std::byte{0x22}, std::byte{0x23}, std::byte{0x24},
                  std::byte{0x25}, std::byte{0x26}, std::byte{0x27}, std::byte{0x28}};
    block.RawSize = static_cast<std::int32_t>(std::size(block.Data));

    CompressCramBlock(block, CramBlockMethod::LZMA);
    EXPECT_EQ(block.Method, CramBlockMethod::LZMA);

    DecompressCramBlock(block);
    EXPECT_EQ(block.Method, CramBlockMethod::RAW);
    EXPECT_EQ(std::size(block.Data), 8u);
    EXPECT_EQ(block.Data[0], std::byte{0x21});
    EXPECT_EQ(block.Data[7], std::byte{0x28});
}
#endif

// ===========================================================================
// Writer / Reader Integration Tests
// ===========================================================================

class CramWriterReaderTest : public ::testing::Test
{
protected:
    tests::TempDirGuard tempDir_;
    std::filesystem::path tmpPath;
    std::filesystem::path tmpCraiPath;

    void SetUp() override
    {
        tempDir_.Reset("cram_writer_reader");
        tmpPath = tempDir_.File("roundtrip.cram");
        tmpCraiPath = DefaultCraiPathFor(tmpPath);
    }

    static SamHeader MakeMinimalHeader()
    {
        SamHeader header;
        header.SetVersion("1.6");
        header.SetSortOrder("unknown");
        header.AddReferenceSequence(ReferenceSequence{"ref", 1000});
        return header;
    }

    static BamRecord MakeUnmappedRecord(std::string_view name, std::string_view seq)
    {
        BamRecord record;
        record.Name(std::string{name})
            .Flag(4)
            .RefId(-1)
            .Pos(-1)
            .MapQ(0)
            .Cigar({})
            .NextRefId(-1)
            .NextPos(-1)
            .Tlen(0)
            .Sequence(std::string{seq});

        std::vector<std::uint8_t> quals(std::size(seq), 30);
        record.Qualities(std::move(quals));
        return record;
    }

    static BamRecord MakeMappedRecord(std::string_view name, std::string_view seq)
    {
        BamRecord record;
        record.Name(std::string{name})
            .Flag(0)
            .RefId(0)
            .Pos(100)
            .MapQ(30)
            .Cigar({CigarOp{CigarOpType::M, static_cast<std::uint32_t>(std::size(seq))}})
            .NextRefId(-1)
            .NextPos(-1)
            .Tlen(0)
            .Sequence(std::string{seq});

        std::vector<std::uint8_t> quals(std::size(seq), 30);
        record.Qualities(std::move(quals));
        return record;
    }
};

TEST(CramWriterConfig, NewFieldsHaveCorrectDefaults)
{
    CramWriterConfig config;
    EXPECT_EQ(config.SlicesPerContainer, 1);
    EXPECT_FALSE(config.CompressionLevel);
    EXPECT_FALSE(config.UseTempFile);
}

TEST_F(CramWriterReaderTest, WriteHeaderOnly)
{
    const SamHeader header = MakeMinimalHeader();
    {
        CramWriter writer{tmpPath, header};
        writer.Close();
    }

    CramReader reader{tmpPath};
    EXPECT_EQ(reader.Header().Version(), "1.6");
    EXPECT_EQ(std::size(reader.Header().ReferenceSequences()), 1u);
    EXPECT_FALSE(reader.ReadRecord());
}

TEST_F(CramWriterReaderTest, WriteAndReadUnmappedRecords)
{
    const SamHeader header = MakeMinimalHeader();
    {
        CramWriter writer{tmpPath, header};
        writer.Write(MakeUnmappedRecord("read1", "ACGT"));
        writer.Write(MakeUnmappedRecord("read2", "TGCA"));
    }

    CramReader reader{tmpPath};
    EXPECT_EQ(reader.Header().Version(), "1.6");

    std::vector<std::string> names;
    std::vector<std::string> seqs;
    for (const auto& record : reader.Records()) {
        names.emplace_back(record.Name());
        seqs.emplace_back(record.Sequence());
    }

    ASSERT_EQ(std::size(names), 2u);
    EXPECT_EQ(names[0], "read1");
    EXPECT_EQ(names[1], "read2");
    EXPECT_EQ(seqs[0], "ACGT");
    EXPECT_EQ(seqs[1], "TGCA");
}

TEST_F(CramWriterReaderTest, ReadRawRecordRoundTrip)
{
    const SamHeader header = MakeMinimalHeader();

    BamRecord source = MakeMappedRecord("raw_decode", "ACGTA");
    source.MutableTags().Set(TagKey{'N', 'M'}, TagValue{std::int64_t{2}});
    {
        CramWriter writer{tmpPath, header};
        writer.Write(source);
    }

    CramReader reader{tmpPath};
    const auto raw = reader.ReadRawRecord();
    ASSERT_TRUE(raw);
    EXPECT_EQ(raw->Name(), "raw_decode");
    EXPECT_EQ(raw->RefId(), 0);
    EXPECT_EQ(raw->Pos(), 100);
    EXPECT_EQ(raw->MapQ(), 30u);
    EXPECT_EQ(raw->Seq().ToString(), "ACGTA");

    const TagMap tags = raw->ParseTags();
    const auto* nm = tags.Get(TagKey{'N', 'M'});
    ASSERT_TRUE(nm);
    const auto* nmValue = std::get_if<std::int64_t>(nm);
    ASSERT_TRUE(nmValue);
    EXPECT_EQ(*nmValue, 2);

    EXPECT_FALSE(reader.ReadRawRecord());
}

TEST_F(CramWriterReaderTest, RawRecordsRangeRoundTrip)
{
    const SamHeader header = MakeMinimalHeader();
    {
        CramWriter writer{tmpPath, header};
        writer.Write(MakeUnmappedRecord("raw_range_1", "AAAA"));
        writer.Write(MakeUnmappedRecord("raw_range_2", "CCCC"));
        writer.Write(MakeUnmappedRecord("raw_range_3", "GGGG"));
    }

    CramReader reader{tmpPath};
    std::vector<std::string> names;
    std::vector<std::string> seqs;
    for (const auto& record : reader.RawRecords()) {
        names.emplace_back(record.Name());
        seqs.emplace_back(record.Seq().ToString());
    }

    ASSERT_EQ(std::size(names), 3u);
    EXPECT_EQ(names[0], "raw_range_1");
    EXPECT_EQ(names[1], "raw_range_2");
    EXPECT_EQ(names[2], "raw_range_3");
    EXPECT_EQ(seqs[0], "AAAA");
    EXPECT_EQ(seqs[1], "CCCC");
    EXPECT_EQ(seqs[2], "GGGG");
}

TEST_F(CramWriterReaderTest, WriteRawRecordRoundTrip)
{
    const SamHeader header = MakeMinimalHeader();

    auto source = MakeMappedRecord("raw_mapped", "ACGTA");
    source.MutableTags().Set(TagKey{'N', 'M'}, TagValue{std::int64_t{1}});
    const std::vector<std::byte> sourceBytes = source.SerializeToBam();
    const RawRecord rawRecord{std::span<const std::byte>{sourceBytes}};

    {
        CramWriter writer{tmpPath, header};
        writer.Write(rawRecord);
    }

    CramReader reader{tmpPath};
    const auto decoded = reader.ReadRecord();
    ASSERT_TRUE(decoded);
    EXPECT_EQ(decoded->Name(), "raw_mapped");
    EXPECT_EQ(decoded->RefId(), 0);
    EXPECT_EQ(decoded->Pos(), 100);
    EXPECT_EQ(decoded->MapQ(), 30u);
    EXPECT_EQ(decoded->Sequence(), "ACGTA");

    const auto* nm = decoded->Tags().Get(TagKey{'N', 'M'});
    ASSERT_TRUE(nm);
    const auto* nmValue = std::get_if<std::int64_t>(nm);
    ASSERT_TRUE(nmValue);
    EXPECT_EQ(*nmValue, 1);
}

TEST_F(CramWriterReaderTest, QueryRawReturnsOverlappingRecords)
{
    const SamHeader header = MakeMinimalHeader();

    CramWriterConfig config;
    config.WriteCrai = true;
    config.RecordsPerSlice = 1;

    {
        CramWriter writer{tmpPath, header, config};

        BamRecord hit = MakeMappedRecord("query_hit", "ACGT");
        hit.Pos(100);
        writer.Write(hit);

        BamRecord miss = MakeMappedRecord("query_miss", "TGCA");
        miss.Pos(180);
        writer.Write(miss);
    }

    const CraiIndex index = CraiIndex::FromFile(tmpCraiPath);
    CramReader reader{tmpPath};
    const auto queried = reader.QueryRaw(index, 0, 100, 110);

    ASSERT_EQ(std::size(queried), 1u);
    EXPECT_EQ(queried[0].Name(), "query_hit");
    EXPECT_EQ(queried[0].Pos(), 100);
    EXPECT_EQ(queried[0].Seq().ToString(), "ACGT");
}

TEST_F(CramWriterReaderTest, WriteRawRecordBatchRoundTrip)
{
    const SamHeader header = MakeMinimalHeader();

    const BamRecord first = MakeUnmappedRecord("batch_read_1", "ACGT");
    BamRecord second = MakeMappedRecord("batch_read_2", "TGCA");
    second.Pos(150);
    second.MutableTags().Set(TagKey{'X', 'Y'}, TagValue{std::string{"value"}});

    std::vector<std::byte> buffer;
    std::vector<RawRecordBatch::RecordExtent> extents;
    AppendSerializedRecord(buffer, extents, first);
    AppendSerializedRecord(buffer, extents, second);

    const RawRecordBatch batch{std::move(buffer), std::move(extents)};

    {
        CramWriter writer{tmpPath, header};
        writer.WriteBatch(batch);
    }

    CramReader reader{tmpPath};
    const auto d1 = reader.ReadRecord();
    const auto d2 = reader.ReadRecord();
    ASSERT_TRUE(d1);
    ASSERT_TRUE(d2);
    EXPECT_EQ(d1->Name(), "batch_read_1");
    EXPECT_EQ(d1->Sequence(), "ACGT");
    EXPECT_EQ(d2->Name(), "batch_read_2");
    EXPECT_EQ(d2->Pos(), 150);
    EXPECT_EQ(d2->Sequence(), "TGCA");

    const auto* xy = d2->Tags().Get(TagKey{'X', 'Y'});
    ASSERT_TRUE(xy);
    const auto* xyValue = std::get_if<std::string>(xy);
    ASSERT_TRUE(xyValue);
    EXPECT_EQ(*xyValue, "value");
}

TEST_F(CramWriterReaderTest, WriteMappedRecordRoundTrip)
{
    const SamHeader header = MakeMinimalHeader();
    {
        CramWriter writer{tmpPath, header};
        writer.Write(MakeMappedRecord("mapped1", "ACGTACGTAC"));
    }

    CramReader reader{tmpPath};
    const auto record = reader.ReadRecord();
    ASSERT_TRUE(record);
    EXPECT_EQ(record->Name(), "mapped1");
    EXPECT_EQ(record->Flag(), 0u);
    EXPECT_EQ(record->RefId(), 0);
    EXPECT_EQ(record->Pos(), 100);
    EXPECT_EQ(record->MapQ(), 30u);
    EXPECT_EQ(record->Sequence(), "ACGTACGTAC");
}

TEST_F(CramWriterReaderTest, CigarRoundTripWithFeatures)
{
    const SamHeader header = MakeMinimalHeader();
    BamRecord record;
    record.Name("cigar_test")
        .Flag(0)
        .RefId(0)
        .Pos(99)
        .MapQ(30)
        .Cigar({
            CigarOp{CigarOpType::S, 3},
            CigarOp{CigarOpType::M, 5},
            CigarOp{CigarOpType::I, 2},
            CigarOp{CigarOpType::M, 3},
            CigarOp{CigarOpType::D, 4},
            CigarOp{CigarOpType::M, 6},
            CigarOp{CigarOpType::S, 3},
        })
        .NextRefId(-1)
        .NextPos(-1)
        .Tlen(0)
        .Sequence("AAACCCCCGGTTTTAAAAAAAC")
        .Qualities(std::vector<std::uint8_t>(22, 30));

    {
        CramWriter writer{tmpPath, header};
        writer.Write(record);
    }

    CramReader reader{tmpPath};
    const auto decoded = reader.ReadRecord();
    ASSERT_TRUE(decoded);

    const auto cigar = decoded->Cigar();
    ASSERT_EQ(std::size(cigar), 7u);
    EXPECT_EQ(cigar[0], CigarOp(CigarOpType::S, 3));
    EXPECT_EQ(cigar[1], CigarOp(CigarOpType::M, 5));
    EXPECT_EQ(cigar[2], CigarOp(CigarOpType::I, 2));
    EXPECT_EQ(cigar[3], CigarOp(CigarOpType::M, 3));
    EXPECT_EQ(cigar[4], CigarOp(CigarOpType::D, 4));
    EXPECT_EQ(cigar[5], CigarOp(CigarOpType::M, 6));
    EXPECT_EQ(cigar[6], CigarOp(CigarOpType::S, 3));
}

TEST_F(CramWriterReaderTest, ReadGroupRoundTrip)
{
    SamHeader header = MakeMinimalHeader();
    header.AddReadGroup(ReadGroup{"rg1"});
    header.AddReadGroup(ReadGroup{"rg2"});

    BamRecord rec1;
    rec1.Name("read1")
        .Flag(0)
        .RefId(0)
        .Pos(100)
        .MapQ(30)
        .Cigar({CigarOp{CigarOpType::M, 5}})
        .NextRefId(-1)
        .NextPos(-1)
        .Tlen(0)
        .Sequence("ACGTG")
        .Qualities({30, 30, 30, 30, 30});
    rec1.MutableTags().Set(TagKey{'R', 'G'}, TagValue{std::string{"rg1"}});

    BamRecord rec2;
    rec2.Name("read2")
        .Flag(0)
        .RefId(0)
        .Pos(200)
        .MapQ(25)
        .Cigar({CigarOp{CigarOpType::M, 4}})
        .NextRefId(-1)
        .NextPos(-1)
        .Tlen(0)
        .Sequence("TTTT")
        .Qualities({20, 20, 20, 20});
    rec2.MutableTags().Set(TagKey{'R', 'G'}, TagValue{std::string{"rg2"}});

    {
        CramWriter writer{tmpPath, header};
        writer.Write(rec1);
        writer.Write(rec2);
    }

    CramReader reader{tmpPath};
    const auto d1 = reader.ReadRecord();
    ASSERT_TRUE(d1);
    const auto* rg1 = d1->Tags().Get(TagKey{'R', 'G'});
    ASSERT_TRUE(rg1);
    const auto* rgText1 = std::get_if<std::string>(rg1);
    ASSERT_TRUE(rgText1);
    EXPECT_EQ(*rgText1, "rg1");

    const auto d2 = reader.ReadRecord();
    ASSERT_TRUE(d2);
    const auto* rg2 = d2->Tags().Get(TagKey{'R', 'G'});
    ASSERT_TRUE(rg2);
    const auto* rgText2 = std::get_if<std::string>(rg2);
    ASSERT_TRUE(rgText2);
    EXPECT_EQ(*rgText2, "rg2");
}

TEST_F(CramWriterReaderTest, TagRoundTrip)
{
    const SamHeader header = MakeMinimalHeader();

    TagArray arrayTag{'C'};
    arrayTag.AppendUInt8(1);
    arrayTag.AppendUInt8(2);
    arrayTag.AppendUInt8(255);

    BamRecord rec;
    rec.Name("tagged")
        .Flag(0)
        .RefId(0)
        .Pos(50)
        .MapQ(30)
        .Cigar({CigarOp{CigarOpType::M, 4}})
        .NextRefId(-1)
        .NextPos(-1)
        .Tlen(0)
        .Sequence("ACGT")
        .Qualities({30, 30, 30, 30});
    rec.MutableTags().Set(TagKey{'T', 'A'}, TagValue{'Q'});
    rec.MutableTags().Set(TagKey{'T', 'c'}, TagValue{std::int64_t{-128}});
    rec.MutableTags().Set(TagKey{'T', 'C'}, TagValue{std::int64_t{255}});
    rec.MutableTags().Set(TagKey{'T', 's'}, TagValue{std::int64_t{-32768}});
    rec.MutableTags().Set(TagKey{'T', 'S'}, TagValue{std::int64_t{65535}});
    rec.MutableTags().Set(TagKey{'T', 'i'}, TagValue{std::int64_t{-2147483648LL}});
    rec.MutableTags().Set(TagKey{'T', 'I'}, TagValue{std::int64_t{4000000000LL}});
    rec.MutableTags().Set(TagKey{'T', 'f'}, TagValue{3.5F});
    rec.MutableTags().Set(TagKey{'T', 'Z'}, TagValue{std::string{"text"}});
    rec.MutableTags().Set(TagKey{'T', 'H'}, TagValue{HexString{"0A0B"}});
    rec.MutableTags().Set(TagKey{'T', 'B'}, TagValue{arrayTag});

    {
        CramWriter writer{tmpPath, header};
        writer.Write(rec);
    }

    CramReader reader{tmpPath};
    const auto decoded = reader.ReadRecord();
    ASSERT_TRUE(decoded);

    const auto* ta = decoded->Tags().Get(TagKey{'T', 'A'});
    ASSERT_TRUE(ta);
    const auto* taChar = std::get_if<char>(ta);
    ASSERT_TRUE(taChar);
    EXPECT_EQ(*taChar, 'Q');

    ExpectDecodedIntTag(*decoded, TagKey{'T', 'c'}, -128);
    ExpectDecodedIntTag(*decoded, TagKey{'T', 'C'}, 255);
    ExpectDecodedIntTag(*decoded, TagKey{'T', 's'}, -32768);
    ExpectDecodedIntTag(*decoded, TagKey{'T', 'S'}, 65535);
    ExpectDecodedIntTag(*decoded, TagKey{'T', 'i'}, -2147483648LL);
    ExpectDecodedIntTag(*decoded, TagKey{'T', 'I'}, 4000000000LL);

    const auto* tf = decoded->Tags().Get(TagKey{'T', 'f'});
    ASSERT_TRUE(tf);
    const auto* tfFloat = std::get_if<float>(tf);
    ASSERT_TRUE(tfFloat);
    EXPECT_FLOAT_EQ(*tfFloat, 3.5F);

    const auto* tz = decoded->Tags().Get(TagKey{'T', 'Z'});
    ASSERT_TRUE(tz);
    const auto* tzText = std::get_if<std::string>(tz);
    ASSERT_TRUE(tzText);
    EXPECT_EQ(*tzText, "text");

    const auto* th = decoded->Tags().Get(TagKey{'T', 'H'});
    ASSERT_TRUE(th);
    const auto* thHex = std::get_if<HexString>(th);
    ASSERT_TRUE(thHex);
    EXPECT_EQ(thHex->value, "0A0B");

    const auto* tb = decoded->Tags().Get(TagKey{'T', 'B'});
    ASSERT_TRUE(tb);
    const auto* tbArray = std::get_if<TagArray>(tb);
    ASSERT_TRUE(tbArray);
    EXPECT_EQ(tbArray->ElementType(), 'C');
    EXPECT_EQ(tbArray->Count(), 3u);
    const auto data = tbArray->Data();
    ASSERT_EQ(std::size(data), 3u);
    EXPECT_EQ(data[0], std::byte{1});
    EXPECT_EQ(data[1], std::byte{2});
    EXPECT_EQ(data[2], std::byte{255});
}

TEST_F(CramWriterReaderTest, WriteMultipleRecords)
{
    const SamHeader header = MakeMinimalHeader();
    {
        CramWriter writer{tmpPath, header};
        for (int i = 0; i < 10; ++i) {
            writer.Write(MakeUnmappedRecord(std::format("read{}", i), "ACGT"));
        }
    }

    CramReader reader{tmpPath};
    std::size_t count = 0;
    for ([[maybe_unused]] const auto& record : reader.Records()) {
        ++count;
    }
    EXPECT_EQ(count, 10u);
}

TEST_F(CramWriterReaderTest, MultiSliceContainerRoundTrip)
{
    const SamHeader header = MakeMinimalHeader();

    CramWriterConfig config;
    config.RecordsPerSlice = 2;
    config.SlicesPerContainer = 3;

    {
        CramWriter writer{tmpPath, header, config};
        for (int i = 0; i < 6; ++i) {
            writer.Write(MakeUnmappedRecord(std::format("read{}", i), "ACGT"));
        }
        writer.Close();
    }

    CramReader reader{tmpPath};
    std::vector<std::string> names;
    for (const auto& record : reader.Records()) {
        names.emplace_back(record.Name());
    }

    ASSERT_EQ(std::size(names), 6u);
    for (int i = 0; i < 6; ++i) {
        EXPECT_EQ(names[i], std::format("read{}", i));
    }
}

TEST_F(CramWriterReaderTest, MultiSliceContainerPartialFlushOnClose)
{
    const SamHeader header = MakeMinimalHeader();

    CramWriterConfig config;
    config.RecordsPerSlice = 3;
    config.SlicesPerContainer = 2;

    {
        CramWriter writer{tmpPath, header, config};
        for (int i = 0; i < 5; ++i) {
            writer.Write(MakeUnmappedRecord(std::format("rec{}", i), "TGCA"));
        }
        writer.Close();
    }

    CramReader reader{tmpPath};
    std::vector<std::string> names;
    for (const auto& record : reader.Records()) {
        names.emplace_back(record.Name());
    }

    ASSERT_EQ(std::size(names), 5u);
}

TEST_F(CramWriterReaderTest, MultiSliceParallelCompressionRoundTrip)
{
    const SamHeader header = MakeMinimalHeader();

    CramWriterConfig config;
    config.RecordsPerSlice = 2;
    config.SlicesPerContainer = 4;
    config.CompressionWorkers = 2;

    {
        CramWriter writer{tmpPath, header, config};
        for (int i = 0; i < 8; ++i) {
            writer.Write(MakeUnmappedRecord(std::format("par{}", i), "ACGTACGT"));
        }
        writer.Close();
    }

    CramReader reader{tmpPath};
    std::vector<std::string> names;
    for (const auto& record : reader.Records()) {
        names.emplace_back(record.Name());
    }

    ASSERT_EQ(std::size(names), 8u);
    for (int i = 0; i < 8; ++i) {
        EXPECT_EQ(names[i], std::format("par{}", i));
    }
}

TEST_F(CramWriterReaderTest, MultiContainerPipelineRoundTrip)
{
    const SamHeader header = MakeMinimalHeader();

    CramWriterConfig config;
    config.RecordsPerSlice = 2;
    config.SlicesPerContainer = 2;
    config.CompressionWorkers = 2;

    {
        CramWriter writer{tmpPath, header, config};
        for (int i = 0; i < 12; ++i) {
            writer.Write(MakeUnmappedRecord(std::format("pipe{}", i), "ACGT"));
        }
        writer.Close();
    }

    CramReader reader{tmpPath};
    std::vector<std::string> names;
    for (const auto& record : reader.Records()) {
        names.emplace_back(record.Name());
    }

    ASSERT_EQ(std::size(names), 12u);
    for (int i = 0; i < 12; ++i) {
        EXPECT_EQ(names[i], std::format("pipe{}", i));
    }
}

TEST_F(CramWriterReaderTest, MultiSliceCraiHasPerSliceEntries)
{
    const SamHeader header = MakeMinimalHeader();

    CramWriterConfig config;
    config.RecordsPerSlice = 1;
    config.SlicesPerContainer = 3;
    config.WriteCrai = true;

    BamRecord r1 = MakeMappedRecord("a", "AAAA");
    r1.Pos(10);
    BamRecord r2 = MakeMappedRecord("b", "CCCC");
    r2.Pos(20);
    BamRecord r3 = MakeMappedRecord("c", "GGGG");
    r3.Pos(30);

    {
        CramWriter writer{tmpPath, header, config};
        writer.Write(r1);
        writer.Write(r2);
        writer.Write(r3);
        writer.Close();
    }

    const CraiIndex index = CraiIndex::FromFile(tmpCraiPath);
    const auto entries = index.EntriesForReference(0);
    ASSERT_EQ(std::size(entries), 3u);

    EXPECT_EQ(entries[0].ContainerOffset, entries[1].ContainerOffset);
    EXPECT_EQ(entries[1].ContainerOffset, entries[2].ContainerOffset);

    EXPECT_LT(entries[0].SliceOffset, entries[1].SliceOffset);
    EXPECT_LT(entries[1].SliceOffset, entries[2].SliceOffset);
}

TEST_F(CramWriterReaderTest, UnlimitedSlicesPerContainerFlushesOnClose)
{
    const SamHeader header = MakeMinimalHeader();

    CramWriterConfig config;
    config.RecordsPerSlice = 1;
    config.SlicesPerContainer = 0;
    config.WriteCrai = true;

    {
        CramWriter writer{tmpPath, header, config};
        writer.Write(MakeMappedRecord("a", "AAAA"));
        writer.Write(MakeMappedRecord("b", "CCCC"));
        writer.Write(MakeMappedRecord("c", "GGGG"));
        writer.Close();
    }

    const CraiIndex index = CraiIndex::FromFile(tmpCraiPath);
    ASSERT_EQ(std::size(index.EntriesForReference(0)), 3u);
    const auto entries = index.EntriesForReference(0);
    EXPECT_EQ(entries[0].ContainerOffset, entries[1].ContainerOffset);
    EXPECT_EQ(entries[1].ContainerOffset, entries[2].ContainerOffset);
}

TEST_F(CramWriterReaderTest, MultiSliceCraiRegionQuery)
{
    const SamHeader header = MakeMinimalHeader();

    CramWriterConfig config;
    config.RecordsPerSlice = 1;
    config.SlicesPerContainer = 3;
    config.WriteCrai = true;

    BamRecord r1 = MakeMappedRecord("early", "AAAAA");
    r1.Pos(10);
    BamRecord r2 = MakeMappedRecord("target", "CCCCC");
    r2.Pos(50);
    BamRecord r3 = MakeMappedRecord("late", "GGGGG");
    r3.Pos(200);

    {
        CramWriter writer{tmpPath, header, config};
        writer.Write(r1);
        writer.Write(r2);
        writer.Write(r3);
        writer.Close();
    }

    const CraiIndex index = CraiIndex::FromFile(tmpCraiPath);
    CramReader reader{tmpPath};
    const auto queried = reader.Query(index, 0, 48, 56);

    ASSERT_EQ(std::size(queried), 1u);
    EXPECT_EQ(queried.front().Name(), "target");
}

TEST_F(CramWriterReaderTest, CompressionLevelRoundTrip)
{
    const SamHeader header = MakeMinimalHeader();

    CramWriterConfig config;
    config.CompressionLevel = 1;

    {
        CramWriter writer{tmpPath, header, config};
        writer.Write(MakeUnmappedRecord("level_test", "ACGTACGTACGT"));
        writer.Close();
    }

    CramReader reader{tmpPath};
    const auto record = reader.ReadRecord();
    ASSERT_TRUE(record);
    EXPECT_EQ(record->Name(), "level_test");
    EXPECT_EQ(record->Sequence(), "ACGTACGTACGT");
}

TEST_F(CramWriterReaderTest, UseTempFileAtomicWrite)
{
    const SamHeader header = MakeMinimalHeader();

    CramWriterConfig config;
    config.UseTempFile = true;

    {
        CramWriter writer{tmpPath, header, config};
        writer.Write(MakeUnmappedRecord("atomic", "ACGT"));
        EXPECT_FALSE(std::filesystem::exists(tmpPath));
        writer.Close();
    }

    EXPECT_TRUE(std::filesystem::exists(tmpPath));

    CramReader reader{tmpPath};
    const auto record = reader.ReadRecord();
    ASSERT_TRUE(record);
    EXPECT_EQ(record->Name(), "atomic");
}

TEST_F(CramWriterReaderTest, GzipCompressionRoundTrip)
{
    const SamHeader header = MakeMinimalHeader();
    {
        CramWriterConfig config;
        config.BlockCompressionMethod = CramBlockMethod::GZIP;
        CramWriter writer{tmpPath, header, config};
        writer.Write(MakeUnmappedRecord("read1", "ACGTACGTACGT"));
    }

    CramReader reader{tmpPath};
    const auto record = reader.ReadRecord();
    ASSERT_TRUE(record);
    EXPECT_EQ(record->Name(), "read1");
    EXPECT_EQ(record->Sequence(), "ACGTACGTACGT");
}

TEST_F(CramWriterReaderTest, RawCompressionRoundTrip)
{
    const SamHeader header = MakeMinimalHeader();
    {
        CramWriterConfig config;
        config.BlockCompressionMethod = CramBlockMethod::RAW;
        CramWriter writer{tmpPath, header, config};
        writer.Write(MakeUnmappedRecord("rawread", "AAAA"));
    }

    CramReader reader{tmpPath};
    const auto record = reader.ReadRecord();
    ASSERT_TRUE(record);
    EXPECT_EQ(record->Name(), "rawread");
    EXPECT_EQ(record->Sequence(), "AAAA");
}

TEST_F(CramWriterReaderTest, UnknownCompressionMethodThrows)
{
    const SamHeader header = MakeMinimalHeader();
    CramWriterConfig config;
    config.BlockCompressionMethod = static_cast<CramBlockMethod>(201);

    EXPECT_THROW(
        WriteSingleRecord(tmpPath, header, config, MakeMappedRecord("codec_record", "ACGTACGT")),
        std::runtime_error);
}

TEST_F(CramWriterReaderTest, Rans4x8CompressionRoundTrip)
{
    const SamHeader header = MakeMinimalHeader();
    {
        CramWriterConfig config;
        config.BlockCompressionMethod = CramBlockMethod::RANS4X8;
        CramWriter writer{tmpPath, header, config};
        writer.Write(MakeUnmappedRecord("r8", "ACGTACGT"));
    }

    CramReader reader{tmpPath};
    const auto record = reader.ReadRecord();
    ASSERT_TRUE(record);
    EXPECT_EQ(record->Name(), "r8");
    EXPECT_EQ(record->Sequence(), "ACGTACGT");
}

TEST_F(CramWriterReaderTest, Rans4x16CompressionRoundTrip)
{
    const SamHeader header = MakeMinimalHeader();
    {
        CramWriterConfig config;
        config.BlockCompressionMethod = CramBlockMethod::RANS4X16;
        CramWriter writer{tmpPath, header, config};
        writer.Write(MakeUnmappedRecord("r16", "TGCATGCA"));
    }

    CramReader reader{tmpPath};
    const auto record = reader.ReadRecord();
    ASSERT_TRUE(record);
    EXPECT_EQ(record->Name(), "r16");
    EXPECT_EQ(record->Sequence(), "TGCATGCA");
}

TEST_F(CramWriterReaderTest, AdaptiveArithCompressionRoundTrip)
{
    const SamHeader header = MakeMinimalHeader();
    {
        CramWriterConfig config;
        config.BlockCompressionMethod = CramBlockMethod::ADAPTIVE_ARITH;
        CramWriter writer{tmpPath, header, config};
        writer.Write(MakeUnmappedRecord("arith", "GATTACA"));
    }

    CramReader reader{tmpPath};
    const auto record = reader.ReadRecord();
    ASSERT_TRUE(record);
    EXPECT_EQ(record->Name(), "arith");
    EXPECT_EQ(record->Sequence(), "GATTACA");
}

TEST_F(CramWriterReaderTest, FqzcompCompressionRoundTrip)
{
    const SamHeader header = MakeMinimalHeader();
    {
        CramWriterConfig config;
        config.BlockCompressionMethod = CramBlockMethod::FQZCOMP;
        CramWriter writer{tmpPath, header, config};
        auto rec = MakeUnmappedRecord("fqz", "ACGTACGT");
        rec.Qualities({10, 11, 12, 13, 14, 15, 16, 17});
        writer.Write(rec);
    }

    CramReader reader{tmpPath};
    const auto record = reader.ReadRecord();
    ASSERT_TRUE(record);
    EXPECT_EQ(record->Name(), "fqz");
    EXPECT_EQ(record->Sequence(), "ACGTACGT");
    ASSERT_EQ(record->Qualities().size(), 8u);
    EXPECT_EQ(record->Qualities()[0], 10);
    EXPECT_EQ(record->Qualities()[7], 17);
}

TEST_F(CramWriterReaderTest, FqzcompRecordsPerSliceRoundTripPreservesQualities)
{
    const SamHeader header = MakeMinimalHeader();
    constexpr std::array<std::int32_t, 3> recordsPerSliceCases{1, 10, 1000};
    constexpr std::array<std::uint16_t, 4> flagCycle{0x41, 0x51, 0x81, 0x91};

    for (const std::int32_t recordsPerSlice : recordsPerSliceCases) {
        SCOPED_TRACE(std::format("RecordsPerSlice={}", recordsPerSlice));

        CramWriterConfig config;
        config.BlockCompressionMethod = CramBlockMethod::FQZCOMP;
        config.RecordsPerSlice = recordsPerSlice;

        const std::size_t totalRecords = FqzcompRoundTripRecordCount(recordsPerSlice);

        std::vector<BamRecord> expectedRecords;
        expectedRecords.reserve(totalRecords);
        {
            CramWriter writer{tmpPath, header, config};
            for (std::size_t i = 0; i < totalRecords; ++i) {
                const auto readLength = 8U + (i % 11);
                std::string sequence;
                sequence.reserve(readLength);
                for (std::size_t j = 0; j < readLength; ++j) {
                    constexpr std::array<char, 4> BASES{'A', 'C', 'G', 'T'};
                    sequence.push_back(BASES[(i + j) % 4]);
                }

                std::vector<std::uint8_t> qualities(readLength);
                for (std::size_t q = 0; q < readLength; ++q) {
                    qualities[q] = static_cast<std::uint8_t>(10 + ((i * 13 + q * 7) % 45));
                }

                auto rec = MakeMappedRecord(std::format("fqz_rec_{:04}", i), sequence);
                rec.Flag(flagCycle[i % (std::size(flagCycle))]);
                rec.Pos(100 + static_cast<std::int32_t>(i * 3));
                rec.Qualities(qualities);

                writer.Write(rec);
                expectedRecords.push_back(std::move(rec));
            }
        }

        CramReader reader{tmpPath};
        std::size_t idx = 0;
        for (const auto& decoded : reader.Records()) {
            ASSERT_LT(idx, std::size(expectedRecords));
            EXPECT_EQ(decoded.Name(), expectedRecords[idx].Name());
            EXPECT_EQ(decoded.Flag(), expectedRecords[idx].Flag());
            const auto decodedQualities = decoded.Qualities();
            const auto expectedQualities = expectedRecords[idx].Qualities();
            ASSERT_EQ(std::size(decodedQualities), std::size(expectedQualities));
            for (std::size_t qi = 0; qi < std::size(decodedQualities); ++qi) {
                EXPECT_EQ(decodedQualities[qi], expectedQualities[qi]);
            }
            ++idx;
        }
        EXPECT_EQ(idx, std::size(expectedRecords));
    }
}

TEST_F(CramWriterReaderTest, NameTokeniserCompressionRoundTrip)
{
    const SamHeader header = MakeMinimalHeader();
    {
        CramWriterConfig config;
        config.BlockCompressionMethod = CramBlockMethod::NAME_TOKENISER;
        CramWriter writer{tmpPath, header, config};
        writer.Write(MakeUnmappedRecord("movie/1/ccs", "ACGT"));
        writer.Write(MakeUnmappedRecord("movie/2/ccs", "TGCA"));
    }

    CramReader reader{tmpPath};
    const auto r1 = reader.ReadRecord();
    const auto r2 = reader.ReadRecord();
    ASSERT_TRUE(r1);
    ASSERT_TRUE(r2);
    EXPECT_EQ(r1->Name(), "movie/1/ccs");
    EXPECT_EQ(r2->Name(), "movie/2/ccs");
}

TEST_F(CramWriterReaderTest, PerDataSeriesCompressionOverridesRoundTrip)
{
    const SamHeader header = MakeMinimalHeader();
    {
        CramWriterConfig config;
        config.BlockCompressionMethod = CramBlockMethod::RAW;
        config.DataSeriesCompressionMethods[CramDataSeries::RN] = CramBlockMethod::NAME_TOKENISER;
        config.DataSeriesCompressionMethods[CramDataSeries::QS] = CramBlockMethod::FQZCOMP;

        CramWriter writer{tmpPath, header, config};

        auto r1 = MakeUnmappedRecord("movie/100/ccs", "ACGTACGT");
        r1.Qualities({10, 11, 12, 13, 14, 15, 16, 17});
        writer.Write(r1);

        auto r2 = MakeUnmappedRecord("movie/101/ccs", "TGCATGCA");
        r2.Qualities({17, 16, 15, 14, 13, 12, 11, 10});
        writer.Write(r2);
    }

    CramReader reader{tmpPath};
    const auto rec1 = reader.ReadRecord();
    const auto rec2 = reader.ReadRecord();
    ASSERT_TRUE(rec1);
    ASSERT_TRUE(rec2);
    EXPECT_EQ(rec1->Name(), "movie/100/ccs");
    EXPECT_EQ(rec2->Name(), "movie/101/ccs");
    ASSERT_EQ(rec1->Qualities().size(), 8u);
    ASSERT_EQ(rec2->Qualities().size(), 8u);
    EXPECT_EQ(rec1->Qualities()[0], 10);
    EXPECT_EQ(rec2->Qualities()[7], 10);
}

TEST_F(CramWriterReaderTest, InvalidPerDataSeriesCompressionOverrideThrows)
{
    const SamHeader header = MakeMinimalHeader();
    CramWriterConfig config;
    config.BlockCompressionMethod = CramBlockMethod::GZIP;
    config.DataSeriesCompressionMethods[CramDataSeries::BA] = CramBlockMethod::FQZCOMP;

    EXPECT_THROW(
        WriteSingleRecord(tmpPath, header, config, MakeUnmappedRecord("invalid_override", "ACGT")),
        std::runtime_error);
}

TEST_F(CramWriterReaderTest, VersionUpgradesTo31ForV31Codecs)
{
    const SamHeader header = MakeMinimalHeader();
    {
        CramWriterConfig config;
        config.BlockCompressionMethod = CramBlockMethod::ADAPTIVE_ARITH;
        config.MinorVersion = 0;
        CramWriter writer{tmpPath, header, config};
        writer.Write(MakeUnmappedRecord("ver", "ACGT"));
    }

    std::ifstream in{tmpPath, std::ios::binary};
    ASSERT_TRUE(in.good());
    std::vector<std::byte> fileDef(26);
    in.read(reinterpret_cast<char*>(fileDef.data()), static_cast<std::streamsize>(fileDef.size()));
    ASSERT_EQ(in.gcount(), static_cast<std::streamsize>(fileDef.size()));

    const auto parsed = ParseFileDefinition(fileDef);
    EXPECT_EQ(parsed.MajorVersion, 3u);
    EXPECT_EQ(parsed.MinorVersion, 1u);
}

TEST_F(CramWriterReaderTest, VersionUpgradesTo31ForPerDataSeriesV31Codec)
{
    const SamHeader header = MakeMinimalHeader();
    {
        CramWriterConfig config;
        config.BlockCompressionMethod = CramBlockMethod::GZIP;
        config.DataSeriesCompressionMethods[CramDataSeries::RN] = CramBlockMethod::NAME_TOKENISER;
        config.MinorVersion = 0;
        CramWriter writer{tmpPath, header, config};
        writer.Write(MakeUnmappedRecord("ver_override", "ACGT"));
    }

    std::ifstream in{tmpPath, std::ios::binary};
    ASSERT_TRUE(in.good());
    std::vector<std::byte> fileDef(26);
    in.read(reinterpret_cast<char*>(fileDef.data()), static_cast<std::streamsize>(fileDef.size()));
    ASSERT_EQ(in.gcount(), static_cast<std::streamsize>(fileDef.size()));

    const auto parsed = ParseFileDefinition(fileDef);
    EXPECT_EQ(parsed.MajorVersion, 3u);
    EXPECT_EQ(parsed.MinorVersion, 1u);
}

TEST_F(CramWriterReaderTest, VersionRemains30WhenGlobalV31MethodIsOverriddenAway)
{
    const SamHeader header = MakeMinimalHeader();
    {
        CramWriterConfig config;
        config.BlockCompressionMethod = CramBlockMethod::NAME_TOKENISER;
        config.DataSeriesCompressionMethods[CramDataSeries::RN] = CramBlockMethod::GZIP;
        config.MinorVersion = 0;
        CramWriter writer{tmpPath, header, config};
        writer.Write(MakeUnmappedRecord("ver_override_down", "ACGT"));
    }

    std::ifstream in{tmpPath, std::ios::binary};
    ASSERT_TRUE(in.good());
    std::vector<std::byte> fileDef(26);
    in.read(reinterpret_cast<char*>(fileDef.data()), static_cast<std::streamsize>(fileDef.size()));
    ASSERT_EQ(in.gcount(), static_cast<std::streamsize>(fileDef.size()));

    const auto parsed = ParseFileDefinition(fileDef);
    EXPECT_EQ(parsed.MajorVersion, 3u);
    EXPECT_EQ(parsed.MinorVersion, 0u);
}

#ifdef PBSAMOA_HAVE_BZIP2
TEST_F(CramWriterReaderTest, Bzip2CompressionRoundTrip)
{
    const SamHeader header = MakeMinimalHeader();
    {
        CramWriterConfig config;
        config.BlockCompressionMethod = CramBlockMethod::BZIP2;
        CramWriter writer{tmpPath, header, config};
        writer.Write(MakeUnmappedRecord("bzread", "ACGTACGT"));
    }

    CramReader reader{tmpPath};
    const auto record = reader.ReadRecord();
    ASSERT_TRUE(record);
    EXPECT_EQ(record->Name(), "bzread");
    EXPECT_EQ(record->Sequence(), "ACGTACGT");
}
#endif

#ifdef PBSAMOA_HAVE_LZMA
TEST_F(CramWriterReaderTest, LzmaCompressionRoundTrip)
{
    const SamHeader header = MakeMinimalHeader();
    {
        CramWriterConfig config;
        config.BlockCompressionMethod = CramBlockMethod::LZMA;
        CramWriter writer{tmpPath, header, config};
        writer.Write(MakeUnmappedRecord("lzread", "TGCATGCA"));
    }

    CramReader reader{tmpPath};
    const auto record = reader.ReadRecord();
    ASSERT_TRUE(record);
    EXPECT_EQ(record->Name(), "lzread");
    EXPECT_EQ(record->Sequence(), "TGCATGCA");
}
#endif

TEST_F(CramWriterReaderTest, QualityScoresPreserved)
{
    const SamHeader header = MakeMinimalHeader();
    {
        CramWriter writer{tmpPath, header};
        auto rec = MakeUnmappedRecord("qread", "ACGT");
        rec.Qualities({10, 20, 30, 40});
        writer.Write(rec);
    }

    CramReader reader{tmpPath};
    const auto record = reader.ReadRecord();
    ASSERT_TRUE(record);
    const auto& quals = record->Qualities();
    ASSERT_EQ(std::size(quals), 4u);
    EXPECT_EQ(quals[0], 10u);
    EXPECT_EQ(quals[1], 20u);
    EXPECT_EQ(quals[2], 30u);
    EXPECT_EQ(quals[3], 40u);
}

TEST_F(CramWriterReaderTest, ThrowOnNonexistent)
{
    EXPECT_THROW(CramReader{"/nonexistent/file.cram"}, std::runtime_error);
}

TEST_F(CramWriterReaderTest, EmptyFileAfterHeader)
{
    const SamHeader header = MakeMinimalHeader();
    {
        const CramWriter writer{tmpPath, header};
        // Write no records
    }

    CramReader reader{tmpPath};
    EXPECT_FALSE(reader.ReadRecord());
}

TEST_F(CramWriterReaderTest, TruncatedContainerThrows)
{
    const SamHeader header = MakeMinimalHeader();
    {
        CramWriter writer{tmpPath, header};
        writer.Write(MakeUnmappedRecord("read1", "ACGT"));
    }

    const auto fileSize = std::filesystem::file_size(tmpPath);
    ASSERT_GT(fileSize, 64u);
    std::filesystem::resize_file(tmpPath, fileSize - 50);  // truncate into final data container

    CramReader reader{tmpPath};
    EXPECT_THROW([[maybe_unused]] auto rec = reader.ReadRecord(), std::runtime_error);
}

TEST_F(CramWriterReaderTest, RangeInterfaceEmptyFile)
{
    const SamHeader header = MakeMinimalHeader();
    {
        const CramWriter writer{tmpPath, header};
    }

    CramReader reader{tmpPath};
    std::size_t count = 0;
    for ([[maybe_unused]] const auto& record : reader.Records()) {
        ++count;
    }
    EXPECT_EQ(count, 0u);
}

TEST_F(CramWriterReaderTest, MoveSemantics)
{
    const SamHeader header = MakeMinimalHeader();
    {
        CramWriter writer{tmpPath, header};
        writer.Write(MakeUnmappedRecord("read1", "ACGT"));

        CramWriter moved{std::move(writer)};
        moved.Write(MakeUnmappedRecord("read2", "TGCA"));
    }

    CramReader reader{tmpPath};
    CramReader moved{std::move(reader)};

    std::size_t count = 0;
    for ([[maybe_unused]] const auto& record : moved.Records()) {
        ++count;
    }
    EXPECT_EQ(count, 2u);
}

TEST_F(CramWriterReaderTest, CraiSidecarDisabledByDefault)
{
    const SamHeader header = MakeMinimalHeader();

    {
        CramWriter writer{tmpPath, header};
        writer.Write(MakeMappedRecord("mapped", "ACGT"));
        writer.Close();
    }

    EXPECT_FALSE(std::filesystem::exists(tmpCraiPath));
}

TEST_F(CramWriterReaderTest, CraiSingleReferenceRowsMatchSliceOffsets)
{
    const SamHeader header = MakeMinimalHeader();

    CramWriterConfig config;
    config.WriteCrai = true;

    {
        CramWriter writer{tmpPath, header, config};
        writer.Write(MakeMappedRecord("mapped1", "ACGTA"));
        writer.Close();
    }

    ASSERT_TRUE(std::filesystem::exists(tmpCraiPath));
    const CraiIndex index = CraiIndex::FromFile(tmpCraiPath);
    ASSERT_EQ(std::size(index.Entries()), 1U);

    const CramSliceOffsetInfo offsets = FirstDataSliceOffsetInfo(tmpPath);
    const CraiEntry& entry = index.Entries().front();
    EXPECT_EQ(entry.SequenceId, 0);
    EXPECT_EQ(entry.AlignmentStart, 101);
    EXPECT_EQ(entry.AlignmentSpan, 5);
    EXPECT_EQ(entry.ContainerOffset, offsets.ContainerOffset);
    EXPECT_EQ(entry.SliceOffset, offsets.SliceOffset);
    EXPECT_EQ(entry.SliceSize, offsets.SliceSize);
}

TEST_F(CramWriterReaderTest, CraiMultiReferenceRowsShareOffsetsAndIncludeUnmapped)
{
    SamHeader header = MakeMinimalHeader();
    header.AddReferenceSequence(ReferenceSequence{"ref2", 1000});

    CramWriterConfig config;
    config.WriteCrai = true;

    BamRecord ref1 = MakeMappedRecord("ref1", "TTTTT");
    ref1.RefId(1).Pos(200);

    BamRecord ref0 = MakeMappedRecord("ref0", "ACGT");
    ref0.RefId(0).Pos(100);

    const BamRecord unmapped = MakeUnmappedRecord("unmapped", "GG");

    {
        CramWriter writer{tmpPath, header, config};
        writer.Write(ref1);
        writer.Write(ref0);
        writer.Write(unmapped);
        writer.Close();
    }

    ASSERT_TRUE(std::filesystem::exists(tmpCraiPath));
    const CraiIndex index = CraiIndex::FromFile(tmpCraiPath);
    ASSERT_EQ(std::size(index.Entries()), 3U);

    const CraiEntry& first = index.Entries().at(0);
    const CraiEntry& second = index.Entries().at(1);
    const CraiEntry& third = index.Entries().at(2);

    EXPECT_EQ(first.SequenceId, 1);
    EXPECT_EQ(first.AlignmentStart, 201);
    EXPECT_EQ(first.AlignmentSpan, 5);

    EXPECT_EQ(second.SequenceId, 0);
    EXPECT_EQ(second.AlignmentStart, 101);
    EXPECT_EQ(second.AlignmentSpan, 4);

    EXPECT_EQ(third.SequenceId, -1);
    EXPECT_EQ(third.AlignmentStart, 0);
    EXPECT_EQ(third.AlignmentSpan, 1);

    EXPECT_EQ(first.ContainerOffset, second.ContainerOffset);
    EXPECT_EQ(first.SliceOffset, second.SliceOffset);
    EXPECT_EQ(first.SliceSize, second.SliceSize);
    EXPECT_EQ(first.ContainerOffset, third.ContainerOffset);
    EXPECT_EQ(first.SliceOffset, third.SliceOffset);
    EXPECT_EQ(first.SliceSize, third.SliceSize);
}

TEST_F(CramWriterReaderTest, CraiRegionQueryMappedReturnsOnlyOverlaps)
{
    const SamHeader header = MakeMinimalHeader();

    CramWriterConfig config;
    config.WriteCrai = true;
    config.RecordsPerSlice = 1;

    BamRecord left = MakeMappedRecord("left", "AAAAA");
    left.Pos(80);
    BamRecord overlap = MakeMappedRecord("overlap", "CCCCC");
    overlap.Pos(100);
    BamRecord right = MakeMappedRecord("right", "GGGGG");
    right.Pos(140);

    {
        CramWriter writer{tmpPath, header, config};
        writer.Write(left);
        writer.Write(overlap);
        writer.Write(right);
        writer.Close();
    }

    const CraiIndex index = CraiIndex::FromFile(tmpCraiPath);
    ASSERT_FALSE(std::empty(index.Entries()));

    std::string duplicatedIndexText;
    for (const CraiEntry& entry : index.Entries()) {
        duplicatedIndexText += std::format(
            "{}\t{}\t{}\t{}\t{}\t{}\n", entry.SequenceId, entry.AlignmentStart, entry.AlignmentSpan,
            entry.ContainerOffset, entry.SliceOffset, entry.SliceSize);
    }
    const CraiEntry& first = index.Entries().front();
    duplicatedIndexText +=
        std::format("{}\t{}\t{}\t{}\t{}\t{}\n", first.SequenceId, first.AlignmentStart,
                    first.AlignmentSpan, first.ContainerOffset, first.SliceOffset, first.SliceSize);

    const std::span<const std::byte> duplicatedTextBytes{
        reinterpret_cast<const std::byte*>(duplicatedIndexText.data()),
        std::size(duplicatedIndexText)};
    const auto duplicatedCompressed = CramGzipCompress(duplicatedTextBytes);
    {
        std::ofstream out{tmpCraiPath, std::ios::binary | std::ios::trunc};
        ASSERT_TRUE(out.is_open());
        out.write(reinterpret_cast<const char*>(duplicatedCompressed.data()),
                  static_cast<std::streamsize>(std::size(duplicatedCompressed)));
        ASSERT_TRUE(out.good());
    }

    const CraiIndex duplicatedIndex = CraiIndex::FromFile(tmpCraiPath);
    CramReader reader{tmpPath};
    const auto queried = reader.Query(duplicatedIndex, 0, 100, 105);

    ASSERT_EQ(std::size(queried), 1U);
    EXPECT_EQ(queried.front().Name(), "overlap");

    for (const BamRecord& record : queried) {
        const std::int32_t recordEnd = QueryRecordEnd(record);
        EXPECT_EQ(record.RefId(), 0);
        EXPECT_LT(record.Pos(), 105);
        EXPECT_GT(recordEnd, 100);
    }
}

TEST_F(CramWriterReaderTest, CraiRegionMultiSliceBoundaryReturnsExpectedReads)
{
    const SamHeader header = MakeMinimalHeader();

    CramWriterConfig config;
    config.WriteCrai = true;
    config.RecordsPerSlice = 1;

    BamRecord leftBoundary = MakeMappedRecord("boundary_left", "AAAAAAA");
    leftBoundary.Pos(95);
    BamRecord rightBoundary = MakeMappedRecord("boundary_right", "CCCCCCC");
    rightBoundary.Pos(104);
    BamRecord outside = MakeMappedRecord("outside", "GGGGGGG");
    outside.Pos(200);

    {
        CramWriter writer{tmpPath, header, config};
        writer.Write(leftBoundary);
        writer.Write(rightBoundary);
        writer.Write(outside);
        writer.Close();
    }

    const CraiIndex index = CraiIndex::FromFile(tmpCraiPath);
    CramReader reader{tmpPath};
    const auto queried = reader.Query(index, 0, 100, 110);

    ASSERT_EQ(std::size(queried), 2U);
    EXPECT_EQ(queried.at(0).Name(), "boundary_left");
    EXPECT_EQ(queried.at(1).Name(), "boundary_right");
}

TEST_F(CramWriterReaderTest, CraiRegionUnmappedQueryReturnsOnlyUnmappedReads)
{
    const SamHeader header = MakeMinimalHeader();

    CramWriterConfig config;
    config.WriteCrai = true;
    config.RecordsPerSlice = 1;

    BamRecord mapped = MakeMappedRecord("mapped", "AAAAA");
    mapped.Pos(100);

    {
        CramWriter writer{tmpPath, header, config};
        writer.Write(mapped);
        writer.Write(MakeUnmappedRecord("unmapped_1", "TTTTT"));
        writer.Write(MakeUnmappedRecord("unmapped_2", "CCCCC"));
        writer.Close();
    }

    const CraiIndex index = CraiIndex::FromFile(tmpCraiPath);
    CramReader reader{tmpPath};
    const auto queried = reader.Query(index, -1, 0, 0);

    ASSERT_EQ(std::size(queried), 2U);
    EXPECT_EQ(queried.at(0).Name(), "unmapped_1");
    EXPECT_EQ(queried.at(1).Name(), "unmapped_2");
    for (const BamRecord& record : queried) {
        EXPECT_LT(record.RefId(), 0);
    }
}

TEST_F(CramWriterReaderTest, CraiRegionEmptyRangeReturnsNoRecords)
{
    const SamHeader header = MakeMinimalHeader();

    CramWriterConfig config;
    config.WriteCrai = true;
    config.RecordsPerSlice = 1;

    BamRecord mapped = MakeMappedRecord("mapped", "AAAAA");
    mapped.Pos(100);

    {
        CramWriter writer{tmpPath, header, config};
        writer.Write(mapped);
        writer.Close();
    }

    const CraiIndex index = CraiIndex::FromFile(tmpCraiPath);
    CramReader reader{tmpPath};
    const auto queried = reader.Query(index, 0, 0, 10);

    EXPECT_TRUE(std::empty(queried));
}

}  // namespace Samoa
}  // namespace PacBio
