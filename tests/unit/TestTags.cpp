#include <pbsamoa/core/Tags.hpp>

#include <gtest/gtest.h>

#include <bit>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {

using ParsedTag = std::optional<std::pair<TagKey, TagValue>>;

TEST(TagKey, ConstructFromChars)
{
    constexpr TagKey key{'N', 'M'};
    EXPECT_EQ(key.First(), 'N');
    EXPECT_EQ(key.Second(), 'M');
}

TEST(TagKey, Equality)
{
    constexpr TagKey a{'N', 'M'};
    constexpr TagKey b{'N', 'M'};
    constexpr TagKey c{'R', 'G'};
    EXPECT_EQ(a, b);
    EXPECT_NE(a, c);
}

TEST(TagKey, ToString)
{
    constexpr TagKey key{'S', 'A'};
    EXPECT_EQ(key.ToString(), "SA");
}

TEST(TagValue, CharType)
{
    const TagValue val{'A'};
    ASSERT_TRUE(std::holds_alternative<char>(val));
    EXPECT_EQ(std::get<char>(val), 'A');
}

TEST(TagValue, IntegerType)
{
    const TagValue val{std::int64_t{42}};
    ASSERT_TRUE(std::holds_alternative<std::int64_t>(val));
    EXPECT_EQ(std::get<std::int64_t>(val), 42);
}

TEST(TagValue, FloatType)
{
    const TagValue val{1.5f};
    ASSERT_TRUE(std::holds_alternative<float>(val));
    EXPECT_FLOAT_EQ(std::get<float>(val), 1.5f);
}

TEST(TagValue, StringType)
{
    const TagValue val{std::string{"hello"}};
    ASSERT_TRUE(std::holds_alternative<std::string>(val));
    EXPECT_EQ(std::get<std::string>(val), "hello");
}

TEST(TagValue, HexStringType)
{
    const TagValue val{HexString{"1AE1"}};
    ASSERT_TRUE(std::holds_alternative<HexString>(val));
    EXPECT_EQ(std::get<HexString>(val).value, "1AE1");
}

TEST(TagArray, Int32Array)
{
    TagArray arr{'i'};
    arr.AppendInt32(10);
    arr.AppendInt32(20);
    arr.AppendInt32(30);

    EXPECT_EQ(arr.Count(), 3U);
    EXPECT_EQ(arr.ElementType(), 'i');
}

TEST(TagMap, SetAndGet)
{
    TagMap tags;
    tags.Set(TagKey{'N', 'M'}, TagValue{std::int64_t{5}});

    const auto* val{tags.Get(TagKey{'N', 'M'})};
    ASSERT_TRUE(val);
    EXPECT_EQ(std::get<std::int64_t>(*val), 5);
}

TEST(TagMap, GetMissing)
{
    const TagMap tags;
    EXPECT_FALSE(tags.Get(TagKey{'N', 'M'}));
}

TEST(TagMap, SetOverwrites)
{
    TagMap tags;
    tags.Set(TagKey{'N', 'M'}, TagValue{std::int64_t{5}});
    tags.Set(TagKey{'N', 'M'}, TagValue{std::int64_t{10}});

    EXPECT_EQ(tags.Size(), 1U);
    EXPECT_EQ(std::get<std::int64_t>(*tags.Get(TagKey{'N', 'M'})), 10);
}

TEST(TagMap, Remove)
{
    TagMap tags;
    tags.Set(TagKey{'N', 'M'}, TagValue{std::int64_t{5}});
    tags.Set(TagKey{'R', 'G'}, TagValue{std::string{"group1"}});

    EXPECT_TRUE(tags.Remove(TagKey{'N', 'M'}));
    EXPECT_EQ(tags.Size(), 1U);
    EXPECT_FALSE(tags.Get(TagKey{'N', 'M'}));

    EXPECT_FALSE(tags.Remove(TagKey{'X', 'X'}));
}

TEST(TagMap, Contains)
{
    TagMap tags;
    tags.Set(TagKey{'N', 'M'}, TagValue{std::int64_t{5}});

    EXPECT_TRUE(tags.Contains(TagKey{'N', 'M'}));
    EXPECT_FALSE(tags.Contains(TagKey{'R', 'G'}));
}

TEST(TagFilter, DropTags)
{
    const DropTags filter{TagKey{'i', 'p'}, TagKey{'p', 'w'}};
    EXPECT_TRUE(filter.ShouldDrop(TagKey{'i', 'p'}));
    EXPECT_TRUE(filter.ShouldDrop(TagKey{'p', 'w'}));
    EXPECT_FALSE(filter.ShouldDrop(TagKey{'N', 'M'}));
}

TEST(TagFilter, KeepTags)
{
    const KeepTags filter{TagKey{'N', 'M'}, TagKey{'R', 'G'}};
    EXPECT_TRUE(filter.ShouldKeep(TagKey{'N', 'M'}));
    EXPECT_TRUE(filter.ShouldKeep(TagKey{'R', 'G'}));
    EXPECT_FALSE(filter.ShouldKeep(TagKey{'i', 'p'}));
}

// --- BAM Tag Parsing Tests ---

namespace {

void AppendU16LE(std::vector<std::byte>& out, std::uint16_t v)
{
    out.push_back(static_cast<std::byte>(v & 0xFFU));
    out.push_back(static_cast<std::byte>((v >> 8U) & 0xFFU));
}

void AppendU32LE(std::vector<std::byte>& out, std::uint32_t v)
{
    out.push_back(static_cast<std::byte>(v & 0xFFU));
    out.push_back(static_cast<std::byte>((v >> 8U) & 0xFFU));
    out.push_back(static_cast<std::byte>((v >> 16U) & 0xFFU));
    out.push_back(static_cast<std::byte>((v >> 24U) & 0xFFU));
}

// Helper: build raw BAM tag bytes for a single integer tag
std::vector<std::byte> MakeBamIntTag(char c1, char c2, char type, std::int64_t value)
{
    std::vector<std::byte> result;
    result.push_back(static_cast<std::byte>(c1));
    result.push_back(static_cast<std::byte>(c2));
    result.push_back(static_cast<std::byte>(type));

    switch (type) {
        case 'c': {
            const std::int8_t v = value;
            result.push_back(static_cast<std::byte>(v));
            break;
        }
        case 'C': {
            const std::uint8_t v = value;
            result.push_back(static_cast<std::byte>(v));
            break;
        }
        case 's': {
            const std::int16_t v{static_cast<std::int16_t>(value)};
            AppendU16LE(result, std::bit_cast<std::uint16_t>(v));
            break;
        }
        case 'S': {
            const std::uint16_t v{static_cast<std::uint16_t>(value)};
            AppendU16LE(result, v);
            break;
        }
        case 'i': {
            const std::int32_t v{static_cast<std::int32_t>(value)};
            AppendU32LE(result, std::bit_cast<std::uint32_t>(v));
            break;
        }
        case 'I': {
            const std::uint32_t v{static_cast<std::uint32_t>(value)};
            AppendU32LE(result, v);
            break;
        }
        default:
            break;
    }
    return result;
}

// Helper: build raw BAM tag bytes for a string tag
std::vector<std::byte> MakeBamStringTag(char c1, char c2, char type, std::string_view value)
{
    std::vector<std::byte> result;
    result.push_back(static_cast<std::byte>(c1));
    result.push_back(static_cast<std::byte>(c2));
    result.push_back(static_cast<std::byte>(type));
    for (const char ch : value) {
        result.push_back(static_cast<std::byte>(ch));
    }
    result.push_back(std::byte{0});  // NUL terminator
    return result;
}

// Helper: build raw BAM tag bytes for a char tag
std::vector<std::byte> MakeBamCharTag(char c1, char c2, char value)
{
    return {static_cast<std::byte>(c1), static_cast<std::byte>(c2), std::byte{'A'},
            static_cast<std::byte>(value)};
}

// Helper: build raw BAM tag bytes for a float tag
std::vector<std::byte> MakeBamFloatTag(char c1, char c2, float value)
{
    std::vector<std::byte> result;
    result.push_back(static_cast<std::byte>(c1));
    result.push_back(static_cast<std::byte>(c2));
    result.push_back(std::byte{'f'});
    AppendU32LE(result, std::bit_cast<std::uint32_t>(value));
    return result;
}

}  // namespace

TEST(TagBamParse, IntegerTypeC)
{
    const std::vector<std::byte> data{MakeBamIntTag('N', 'M', 'C', 5)};
    const TagMap tags{ParseTagsFromBam(data)};
    ASSERT_EQ(tags.Size(), 1U);
    const auto* val{tags.Get(TagKey{'N', 'M'})};
    ASSERT_TRUE(val);
    EXPECT_EQ(std::get<std::int64_t>(*val), 5);
}

TEST(TagBamParse, IntegerTypeI)
{
    const std::vector<std::byte> data{MakeBamIntTag('X', 'Y', 'I', 100000)};
    const TagMap tags{ParseTagsFromBam(data)};
    const auto* val{tags.Get(TagKey{'X', 'Y'})};
    ASSERT_TRUE(val);
    EXPECT_EQ(std::get<std::int64_t>(*val), 100000);
}

TEST(TagBamParse, SignedInt)
{
    const std::vector<std::byte> data{MakeBamIntTag('A', 'S', 'c', -10)};
    const TagMap tags{ParseTagsFromBam(data)};
    const auto* val{tags.Get(TagKey{'A', 'S'})};
    ASSERT_TRUE(val);
    EXPECT_EQ(std::get<std::int64_t>(*val), -10);
}

TEST(TagBamParse, UnknownBArraySubtypeStopsParsing)
{
    // 'XX:B:x' has an unknown element subtype, so its payload length is undefined; parsing
    // must stop instead of mis-syncing onto the following bytes (htslib bad_aux behaviour).
    const std::vector<std::byte> data{
        std::byte{'X'},  std::byte{'X'}, std::byte{'B'}, std::byte{'x'},
        std::byte{1},    std::byte{0},   std::byte{0},   std::byte{0},  // count = 1
        std::byte{0xAB},                                                // bogus payload byte
        std::byte{'N'},  std::byte{'M'}, std::byte{'C'}, std::byte{7},  // must NOT be reached
    };
    const TagMap tags{ParseTagsFromBam(data)};
    EXPECT_FALSE(tags.Contains(TagKey{'X', 'X'}));
    EXPECT_FALSE(tags.Contains(TagKey{'N', 'M'}));
}

TEST(TagBamParse, StringTypeZ)
{
    const std::vector<std::byte> data{MakeBamStringTag('R', 'G', 'Z', "group1")};
    const TagMap tags{ParseTagsFromBam(data)};
    const auto* val{tags.Get(TagKey{'R', 'G'})};
    ASSERT_TRUE(val);
    EXPECT_EQ(std::get<std::string>(*val), "group1");
}

TEST(TagBamParse, CharTypeA)
{
    const std::vector<std::byte> data{MakeBamCharTag('X', 'S', '+')};
    const TagMap tags{ParseTagsFromBam(data)};
    const auto* val{tags.Get(TagKey{'X', 'S'})};
    ASSERT_TRUE(val);
    EXPECT_EQ(std::get<char>(*val), '+');
}

TEST(TagBamParse, FloatType)
{
    const std::vector<std::byte> data{MakeBamFloatTag('Z', 'S', 3.14f)};
    const TagMap tags{ParseTagsFromBam(data)};
    const auto* val{tags.Get(TagKey{'Z', 'S'})};
    ASSERT_TRUE(val);
    EXPECT_FLOAT_EQ(std::get<float>(*val), 3.14f);
}

TEST(TagBamParse, MultipleTags)
{
    std::vector<std::byte> data{MakeBamIntTag('N', 'M', 'C', 5)};
    const std::vector<std::byte> rg{MakeBamStringTag('R', 'G', 'Z', "lane1")};
    data.insert(std::ranges::end(data), std::ranges::begin(rg), std::ranges::end(rg));

    const TagMap tags{ParseTagsFromBam(data)};
    EXPECT_EQ(tags.Size(), 2U);
    EXPECT_EQ(std::get<std::int64_t>(*tags.Get(TagKey{'N', 'M'})), 5);
    EXPECT_EQ(std::get<std::string>(*tags.Get(TagKey{'R', 'G'})), "lane1");
}

TEST(TagBamParse, ArrayTypeBI)
{
    // Build B:I tag with 3 elements: [10, 20, 30]
    std::vector<std::byte> data;
    data.push_back(static_cast<std::byte>('X'));
    data.push_back(static_cast<std::byte>('A'));
    data.push_back(static_cast<std::byte>('B'));
    data.push_back(static_cast<std::byte>('I'));

    const std::uint32_t count{3};
    AppendU32LE(data, count);

    for (const std::uint32_t val : {10U, 20U, 30U}) {
        AppendU32LE(data, val);
    }

    const TagMap tags{ParseTagsFromBam(data)};
    const auto* val{tags.Get(TagKey{'X', 'A'})};
    ASSERT_TRUE(val);
    ASSERT_TRUE(std::holds_alternative<TagArray>(*val));

    const TagArray& arr{std::get<TagArray>(*val)};
    EXPECT_EQ(arr.ElementType(), 'I');
    EXPECT_EQ(arr.Count(), 3U);
}

TEST(TagBamParse, EmptyData)
{
    const TagMap tags{ParseTagsFromBam({})};
    EXPECT_TRUE(tags.Empty());
}

// --- SAM Tag Parsing Tests ---

TEST(TagSamParse, RejectsInvalidTagKey)
{
    // SAMv1 §1.5: TAG must match [A-Za-z][A-Za-z0-9].
    EXPECT_FALSE(ParseTagFromSam("1X:i:5").has_value());  // first char must be a letter
    EXPECT_FALSE(ParseTagFromSam("X-:i:5").has_value());  // second char must be alphanumeric
}

TEST(TagSamParse, RejectsNonPrintableTypeAValue)
{
    // SAMv1 §1.5: type 'A' must be a printable character [!-~]; space (0x20) is excluded.
    EXPECT_FALSE(ParseTagFromSam("XA:A: ").has_value());
    EXPECT_TRUE(ParseTagFromSam("XA:A:Q").has_value());
}

TEST(TagSamParse, RejectsTypeZWithEmbeddedTab)
{
    // SAMv1 §1.5: type 'Z' is [ !-~]*; an embedded tab would corrupt the SAM column layout.
    EXPECT_FALSE(ParseTagFromSam(std::string{"XZ:Z:ab"} + '\t' + "cd").has_value());
    EXPECT_TRUE(ParseTagFromSam("XZ:Z:ab cd").has_value());  // space is permitted
}

TEST(TagBamParse, RejectsMalformedHexHtag)
{
    // SAMv1 §4.2.4: 'H' is an even-length run of hex digits. The BAM decode path must
    // reject malformed hex (matching the SAM-text path), not round-trip garbage.
    // Aux layout: 2-byte tag key, 1-byte type 'H', then NUL-terminated hex text.
    const auto auxBytes{[](std::string_view body) {
        std::vector<std::byte> v;
        for (const char c : body) {
            v.push_back(static_cast<std::byte>(static_cast<unsigned char>(c)));
        }
        v.push_back(std::byte{0});  // NUL terminator
        return v;
    }};
    EXPECT_THROW(ParseTagsFromBam(auxBytes("XHHABC")), std::runtime_error);  // 3 digits: odd
    EXPECT_THROW(ParseTagsFromBam(auxBytes("XHHAZ")), std::runtime_error);   // 'Z' not hex
    EXPECT_NO_THROW(ParseTagsFromBam(auxBytes("XHHAB")));                    // valid: "AB"
}

TEST(TagSamParse, IntegerTag)
{
    const ParsedTag result{ParseTagFromSam("NM:i:5")};
    ASSERT_TRUE(result);
    EXPECT_EQ(result->first, TagKey('N', 'M'));
    EXPECT_EQ(std::get<std::int64_t>(result->second), 5);
}

TEST(TagSamParse, NegativeInteger)
{
    const ParsedTag result{ParseTagFromSam("AS:i:-10")};
    ASSERT_TRUE(result);
    EXPECT_EQ(std::get<std::int64_t>(result->second), -10);
}

TEST(TagSamParse, IntegerOutOfBamRangeRejected)
{
    // [INT32_MIN, UINT32_MAX] is the representable range; htslib rejects wider values.
    EXPECT_TRUE(ParseTagFromSam("ZZ:i:4294967295").has_value());    // UINT32_MAX: ok
    EXPECT_FALSE(ParseTagFromSam("ZZ:i:4294967296").has_value());   // > UINT32_MAX
    EXPECT_TRUE(ParseTagFromSam("ZZ:i:-2147483648").has_value());   // INT32_MIN: ok
    EXPECT_FALSE(ParseTagFromSam("ZZ:i:-2147483649").has_value());  // < INT32_MIN
}

TEST(TagSerialize, IntegerOutOfBamRangeThrows)
{
    // The BAM encode path must fail loud rather than silently truncate to 32 bits.
    TagMap tags;
    tags.Set(TagKey{'Z', 'Z'}, TagValue{std::int64_t{5'000'000'000}});  // > UINT32_MAX
    EXPECT_THROW((void)SerializeTagsToBam(tags), std::runtime_error);
}

TEST(TagSamParse, StringTag)
{
    const ParsedTag result{ParseTagFromSam("RG:Z:lane1")};
    ASSERT_TRUE(result);
    EXPECT_EQ(std::get<std::string>(result->second), "lane1");
}

TEST(TagSamParse, StringWithColons)
{
    const ParsedTag result{ParseTagFromSam("SA:Z:ref,29,-,6H5M,17,0;")};
    ASSERT_TRUE(result);
    EXPECT_EQ(std::get<std::string>(result->second), "ref,29,-,6H5M,17,0;");
}

TEST(TagSamParse, CharTag)
{
    const ParsedTag result{ParseTagFromSam("XS:A:+")};
    ASSERT_TRUE(result);
    EXPECT_EQ(std::get<char>(result->second), '+');
}

TEST(TagSamParse, FloatTag)
{
    const ParsedTag result{ParseTagFromSam("ZS:f:3.14")};
    ASSERT_TRUE(result);
    EXPECT_NEAR(std::get<float>(result->second), 3.14f, 0.001f);
}

TEST(TagSamParse, ArrayTag)
{
    const ParsedTag result{ParseTagFromSam("XA:B:i,10,20,30")};
    ASSERT_TRUE(result);
    const TagArray& arr{std::get<TagArray>(result->second)};
    EXPECT_EQ(arr.ElementType(), 'i');
    EXPECT_EQ(arr.Count(), 3U);
}

TEST(TagSamParse, InvalidFormat)
{
    EXPECT_FALSE(ParseTagFromSam("invalid"));
    EXPECT_FALSE(ParseTagFromSam("NM"));
    EXPECT_FALSE(ParseTagFromSam(""));
}

TEST(TagSamParse, RejectsTrailingJunkInNumericValue)
{
    EXPECT_FALSE(ParseTagFromSam("NM:i:5x"));
    EXPECT_FALSE(ParseTagFromSam("ZS:f:3.14abc"));
}

TEST(TagSamParse, RejectsInvalidArrayElements)
{
    EXPECT_FALSE(ParseTagFromSam("XA:B:i,10,20x,30"));
    EXPECT_FALSE(ParseTagFromSam("XA:B:C,10,999"));
    EXPECT_FALSE(ParseTagFromSam("XA:B:Q,1,2"));
    EXPECT_FALSE(ParseTagFromSam("XA:B:i,"));
}

// --- SAM Tag Serialization Tests ---

TEST(TagSerialize, ToSamInteger)
{
    const std::string result{SerializeTagToSam(TagKey{'N', 'M'}, TagValue{std::int64_t{5}})};
    EXPECT_EQ(result, "NM:i:5");
}

TEST(TagSerialize, ToSamString)
{
    const std::string result{SerializeTagToSam(TagKey{'R', 'G'}, TagValue{std::string{"lane1"}})};
    EXPECT_EQ(result, "RG:Z:lane1");
}

TEST(TagSerialize, ToSamChar)
{
    const std::string result{SerializeTagToSam(TagKey{'X', 'S'}, TagValue{'+'})};
    EXPECT_EQ(result, "XS:A:+");
}

TEST(TagSerialize, BamIntegerUsesLittleEndianEncoding)
{
    TagMap tags;
    tags.Set(TagKey{'N', 'M'}, TagValue{std::int64_t{0x1234}});

    const std::vector<std::byte> bam{SerializeTagsToBam(tags)};
    ASSERT_EQ(std::size(bam), 5U);
    EXPECT_EQ(bam[0], static_cast<std::byte>('N'));
    EXPECT_EQ(bam[1], static_cast<std::byte>('M'));
    EXPECT_EQ(bam[2], static_cast<std::byte>('S'));
    EXPECT_EQ(bam[3], static_cast<std::byte>(0x34));
    EXPECT_EQ(bam[4], static_cast<std::byte>(0x12));
}

TEST(TagSerialize, BamFloatUsesLittleEndianEncoding)
{
    TagMap tags;
    tags.Set(TagKey{'Z', 'S'}, TagValue{1.0F});

    const std::vector<std::byte> bam{SerializeTagsToBam(tags)};
    ASSERT_EQ(std::size(bam), 7U);
    EXPECT_EQ(bam[0], static_cast<std::byte>('Z'));
    EXPECT_EQ(bam[1], static_cast<std::byte>('S'));
    EXPECT_EQ(bam[2], static_cast<std::byte>('f'));
    EXPECT_EQ(bam[3], static_cast<std::byte>(0x00));
    EXPECT_EQ(bam[4], static_cast<std::byte>(0x00));
    EXPECT_EQ(bam[5], static_cast<std::byte>(0x80));
    EXPECT_EQ(bam[6], static_cast<std::byte>(0x3F));
}

// --- BAM Round-Trip Tests ---

TEST(TagSerialize, BamRoundTripInteger)
{
    TagMap original;
    original.Set(TagKey{'N', 'M'}, TagValue{std::int64_t{5}});
    original.Set(TagKey{'A', 'S'}, TagValue{std::int64_t{-42}});
    original.Set(TagKey{'X', 'I'}, TagValue{std::int64_t{100000}});

    const std::vector<std::byte> bamBytes{SerializeTagsToBam(original)};
    const TagMap parsed{ParseTagsFromBam(bamBytes)};

    EXPECT_EQ(parsed.Size(), 3U);
    EXPECT_EQ(std::get<std::int64_t>(*parsed.Get(TagKey{'N', 'M'})), 5);
    EXPECT_EQ(std::get<std::int64_t>(*parsed.Get(TagKey{'A', 'S'})), -42);
    EXPECT_EQ(std::get<std::int64_t>(*parsed.Get(TagKey{'X', 'I'})), 100000);
}

TEST(TagSerialize, BamRoundTripMixed)
{
    TagMap original;
    original.Set(TagKey{'N', 'M'}, TagValue{std::int64_t{1}});
    original.Set(TagKey{'R', 'G'}, TagValue{std::string{"group1"}});
    original.Set(TagKey{'X', 'S'}, TagValue{'+'});
    original.Set(TagKey{'Z', 'S'}, TagValue{1.5f});

    const std::vector<std::byte> bamBytes{SerializeTagsToBam(original)};
    const TagMap parsed{ParseTagsFromBam(bamBytes)};

    EXPECT_EQ(parsed.Size(), 4U);
    EXPECT_EQ(std::get<std::int64_t>(*parsed.Get(TagKey{'N', 'M'})), 1);
    EXPECT_EQ(std::get<std::string>(*parsed.Get(TagKey{'R', 'G'})), "group1");
    EXPECT_EQ(std::get<char>(*parsed.Get(TagKey{'X', 'S'})), '+');
    EXPECT_FLOAT_EQ(std::get<float>(*parsed.Get(TagKey{'Z', 'S'})), 1.5f);
}

TEST(TagSamParse, HexStringTag)
{
    const ParsedTag result{ParseTagFromSam("BC:H:1AE1")};
    ASSERT_TRUE(result);
    EXPECT_EQ(result->first, TagKey('B', 'C'));
    ASSERT_TRUE(std::holds_alternative<HexString>(result->second));
    EXPECT_EQ(std::get<HexString>(result->second).value, "1AE1");
}

TEST(TagSamParse, RejectsInvalidHexStringTag)
{
    EXPECT_FALSE(ParseTagFromSam("BC:H:ABC"));   // odd length
    EXPECT_FALSE(ParseTagFromSam("BC:H:1AG1"));  // non-hex char
}

TEST(TagSerialize, ToSamHexString)
{
    const std::string result{SerializeTagToSam(TagKey{'B', 'C'}, TagValue{HexString{"1AE1"}})};
    EXPECT_EQ(result, "BC:H:1AE1");
}

TEST(TagSerialize, BamRoundTripHexString)
{
    TagMap original;
    original.Set(TagKey{'B', 'C'}, TagValue{HexString{"1AE1"}});

    const std::vector<std::byte> bamBytes{SerializeTagsToBam(original)};
    const TagMap parsed{ParseTagsFromBam(bamBytes)};

    ASSERT_EQ(parsed.Size(), 1U);
    const auto* val{parsed.Get(TagKey{'B', 'C'})};
    ASSERT_TRUE(val);
    ASSERT_TRUE(std::holds_alternative<HexString>(*val));
    EXPECT_EQ(std::get<HexString>(*val).value, "1AE1");
}

}  // namespace Samoa
}  // namespace PacBio
