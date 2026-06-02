#include <pbsamoa/core/LosslessClipping.hpp>

#include <pbsamoa/core/BamRecord.hpp>
#include <pbsamoa/core/Tags.hpp>

#include <gtest/gtest.h>

#include <algorithm>
#include <ranges>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {
namespace {

TagArray MakeUInt8Array(const std::vector<std::uint8_t>& values)
{
    TagArray array{'C'};
    for (const std::uint8_t value : values) {
        array.AppendUInt8(value);
    }
    return array;
}

std::vector<std::uint8_t> Ramp(std::uint8_t start, std::int32_t count)
{
    std::vector<std::uint8_t> values;
    values.reserve(count);
    for (std::int32_t i = 0; i < count; ++i) {
        values.push_back(static_cast<std::uint8_t>(start + i));
    }
    return values;
}

const TagArray& ArrayTag(const BamRecord& record, TagKey key)
{
    return std::get<TagArray>(*record.Tags().Get(key));
}

// Decode a 'C' (uint8) array to readable values for diff-on-failure.
std::vector<std::uint8_t> DecodeC(const TagArray& array)
{
    std::vector<std::uint8_t> out;
    for (const std::byte b : array.Data()) {
        out.push_back(std::to_integer<std::uint8_t>(b));
    }
    return out;
}

std::vector<std::uint8_t> TagValues(const BamRecord& record, TagKey key)
{
    return DecodeC(ArrayTag(record, key));
}

// Unmapped CCS-style record with distinct forward (fp/fi) and reverse (rp/ri)
// per-base tracks so a mishandled orientation would be caught.
BamRecord MakeUnmappedRecord(std::int32_t n)
{
    std::string seq(static_cast<std::size_t>(n), 'A');
    for (std::int32_t i = 0; i < n; ++i) {
        seq[i] = "ACGT"[i % 4];
    }

    BamRecord record;
    record.Name("movie/1/ccs").Flag(0x4).Sequence(seq).Qualities(Ramp(20, n));

    TagMap tags;
    tags.Set(TagKey{'f', 'p'}, MakeUInt8Array(Ramp(0, n)));
    tags.Set(TagKey{'f', 'i'}, MakeUInt8Array(Ramp(40, n)));
    tags.Set(TagKey{'r', 'p'}, MakeUInt8Array(Ramp(100, n)));
    tags.Set(TagKey{'r', 'i'}, MakeUInt8Array(Ramp(140, n)));
    record.Tags(std::move(tags));
    return record;
}

}  // namespace

TEST(LosslessClipping, RoundTripRestoresOriginal)
{
    constexpr std::int32_t n{20};
    const BamRecord original{MakeUnmappedRecord(n)};

    BamRecord record{MakeUnmappedRecord(n)};
    ClipToQueryLossless(record, 4, 16);

    // Clipped to the insert, with a restore blob attached.
    EXPECT_EQ(std::ssize(record.Sequence()), 12);
    EXPECT_EQ(ArrayTag(record, TagKey{'f', 'p'}).Count(), 12U);
    EXPECT_EQ(ArrayTag(record, TagKey{'r', 'p'}).Count(), 12U);
    EXPECT_TRUE(record.Tags().Contains(TagKey{'l', 's'}));

    ASSERT_TRUE(RestoreFromLossless(record));

    EXPECT_EQ(record.Sequence(), original.Sequence());
    ASSERT_EQ(std::ssize(record.Qualities()), n);
    EXPECT_TRUE(std::ranges::equal(record.Qualities(), original.Qualities()));
    EXPECT_EQ(TagValues(record, TagKey{'f', 'p'}), TagValues(original, TagKey{'f', 'p'}));
    EXPECT_EQ(TagValues(record, TagKey{'f', 'i'}), TagValues(original, TagKey{'f', 'i'}));
    EXPECT_EQ(TagValues(record, TagKey{'r', 'p'}), TagValues(original, TagKey{'r', 'p'}));
    EXPECT_EQ(TagValues(record, TagKey{'r', 'i'}), TagValues(original, TagKey{'r', 'i'}));
    EXPECT_FALSE(record.Tags().Contains(TagKey{'l', 's'}));
}

TEST(LosslessClipping, NestedClipsFullyRestore)
{
    constexpr std::int32_t n{20};
    const BamRecord original{MakeUnmappedRecord(n)};

    BamRecord record{MakeUnmappedRecord(n)};
    ClipToQueryLossless(record, 4, 16);  // 20 -> 12, keep [4,16)
    ClipToQueryLossless(record, 2, 10);  // 12 -> 8,  keep [2,10) of the clipped read
    EXPECT_EQ(std::ssize(record.Sequence()), 8);

    ASSERT_TRUE(RestoreFromLossless(record));  // unwinds both nested levels

    EXPECT_EQ(record.Sequence(), original.Sequence());
    EXPECT_TRUE(std::ranges::equal(record.Qualities(), original.Qualities()));
    EXPECT_EQ(ArrayTag(record, TagKey{'f', 'p'}), ArrayTag(original, TagKey{'f', 'p'}));
    EXPECT_EQ(ArrayTag(record, TagKey{'r', 'p'}), ArrayTag(original, TagKey{'r', 'p'}));
    EXPECT_FALSE(record.Tags().Contains(TagKey{'l', 's'}));
}

// MM/ML survive a round trip when modifications straddle the lead, retained, and
// trailing regions. Sequence "ACGT"*5 has C at positions 1,5,9,13,17; the MM
// "C+m?,0,1,1" marks C-indices 0,2,4 (positions 1,9,17), i.e. one mod in each of
// lead [0,4), retained [4,16), trailing [16,20).
TEST(LosslessClipping, RoundTripRestoresBaseModsAcrossFlanks)
{
    BamRecord record{MakeUnmappedRecord(20)};
    TagMap tags{record.Tags()};
    tags.Set(TagKey{'M', 'M'}, std::string{"C+m?,0,1,1;"});
    tags.Set(TagKey{'M', 'L'}, MakeUInt8Array({200, 150, 100}));
    record.Tags(std::move(tags));

    ClipToQueryLossless(record, 4, 16);
    ASSERT_TRUE(RestoreFromLossless(record));

    const TagValue* mm{record.Tags().Get(TagKey{'M', 'M'})};
    ASSERT_NE(mm, nullptr);
    EXPECT_EQ(std::get<std::string>(*mm), "C+m?,0,1,1;");
    EXPECT_EQ(TagValues(record, TagKey{'M', 'L'}), (std::vector<std::uint8_t>{200, 150, 100}));
}

// All modifications fall inside the retained window (none in the flanks), so the
// first retained skip is rewritten on clip and must be recovered from pMM.
// "C+m?,1,1" marks C-indices 1,3 (positions 5,13), both inside [4,16).
TEST(LosslessClipping, RoundTripRecoversFirstRetainedSkipViaPrefixFix)
{
    BamRecord record{MakeUnmappedRecord(20)};
    TagMap tags{record.Tags()};
    tags.Set(TagKey{'M', 'M'}, std::string{"C+m?,1,1;"});
    tags.Set(TagKey{'M', 'L'}, MakeUInt8Array({90, 80}));
    record.Tags(std::move(tags));

    ClipToQueryLossless(record, 4, 16);
    // The clipped read's first skip was rewritten (a C at position 1 is now gone).
    ASSERT_TRUE(RestoreFromLossless(record));

    const TagValue* mm{record.Tags().Get(TagKey{'M', 'M'})};
    ASSERT_NE(mm, nullptr);
    EXPECT_EQ(std::get<std::string>(*mm), "C+m?,1,1;");
    EXPECT_EQ(TagValues(record, TagKey{'M', 'L'}), (std::vector<std::uint8_t>{90, 80}));
}

// All modifications fall in the flanks with none retained, while canonical bases
// before the clip are dropped (pMM > 0). The lead+trail skips must reconstruct the
// original exactly without the (unused) prefix fixup leaking into the trail.
// Sequence "CCAAA...CC": C at positions 0,1,18,19; "C+m?,0,1" marks C-indices 0,2
// (positions 0 and 18) — lead and trail, nothing in [4,16).
TEST(LosslessClipping, RoundTripBaseModsWithNoneRetained)
{
    std::string seq(20, 'A');
    seq[0] = 'C';
    seq[1] = 'C';
    seq[18] = 'C';
    seq[19] = 'C';

    BamRecord record;
    record.Name("movie/1/ccs").Flag(0x4).Sequence(seq).Qualities(Ramp(20, 20));
    TagMap tags;
    tags.Set(TagKey{'M', 'M'}, std::string{"C+m?,0,1;"});
    tags.Set(TagKey{'M', 'L'}, MakeUInt8Array({200, 100}));
    record.Tags(std::move(tags));

    ClipToQueryLossless(record, 4, 16);
    ASSERT_TRUE(RestoreFromLossless(record));

    const TagValue* mm{record.Tags().Get(TagKey{'M', 'M'})};
    ASSERT_NE(mm, nullptr);
    EXPECT_EQ(std::get<std::string>(*mm), "C+m?,0,1;");
    EXPECT_EQ(TagValues(record, TagKey{'M', 'L'}), (std::vector<std::uint8_t>{200, 100}));
}

// Nested lossless clips with base mods: two successive clips, each splitting the
// MM/ML calls across lead/retained/trail, must both unwind to the original. The
// inner clip's pMM fixup has to compose with the outer clip's on restore.
// Sequence "ACGT"*5 has C at positions 1,5,9,13,17; "C+m?,0,1,1" marks C-indices
// 0,2,4 (positions 1,9,17) — only the position-9 mod survives the first clip.
TEST(LosslessClipping, NestedClipsRestoreBaseMods)
{
    BamRecord record{MakeUnmappedRecord(20)};
    TagMap tags{record.Tags()};
    tags.Set(TagKey{'M', 'M'}, std::string{"C+m?,0,1,1;"});
    tags.Set(TagKey{'M', 'L'}, MakeUInt8Array({200, 150, 100}));
    record.Tags(std::move(tags));

    ClipToQueryLossless(record, 4, 16);  // 20 -> 12, mod at position 9 retained
    ClipToQueryLossless(record, 2, 10);  // 12 -> 8,  mod still retained
    ASSERT_TRUE(RestoreFromLossless(record));

    const TagValue* mm{record.Tags().Get(TagKey{'M', 'M'})};
    ASSERT_NE(mm, nullptr);
    EXPECT_EQ(std::get<std::string>(*mm), "C+m?,0,1,1;");
    EXPECT_EQ(TagValues(record, TagKey{'M', 'L'}), (std::vector<std::uint8_t>{200, 150, 100}));
}

// A multi-code prefix ("C+mh") carries ModCodeCount==2 interleaved ML values per
// site. A site fully inside the retained window must keep BOTH probabilities across
// a clip/restore — the restore-path ML seeding has to advance by skip-count * stride,
// not skip-count. Sequence "ACGT"*5 has C at positions 1,5,9,13,17; "C+mh,2" marks
// C-index 2 (position 9), inside [4,16).
TEST(LosslessClipping, RoundTripRestoresMultiCodeBaseModMlValues)
{
    BamRecord record{MakeUnmappedRecord(20)};
    TagMap tags{record.Tags()};
    tags.Set(TagKey{'M', 'M'}, std::string{"C+mh,2;"});
    tags.Set(TagKey{'M', 'L'}, MakeUInt8Array({210, 110}));  // {m, h} for the single site
    record.Tags(std::move(tags));

    ClipToQueryLossless(record, 4, 16);
    ASSERT_TRUE(RestoreFromLossless(record));

    const TagValue* mm{record.Tags().Get(TagKey{'M', 'M'})};
    ASSERT_NE(mm, nullptr);
    EXPECT_EQ(std::get<std::string>(*mm), "C+mh,2;");
    EXPECT_EQ(TagValues(record, TagKey{'M', 'L'}), (std::vector<std::uint8_t>{210, 110}));
}

// Canonical base 'N' is the MM wildcard: it matches every base, not the literal
// character 'N'. On an A/C/G/T read (no literal N) with an "N+" mod and a left clip,
// the prefix-lost fixup must use the wildcard base count so the first retained skip
// restores to its original value. "N+m?,5,3" marks bases 5 and 9, both inside [4,16).
TEST(LosslessClipping, RoundTripRecoversWildcardPrefixSkip)
{
    BamRecord record{MakeUnmappedRecord(20)};  // "ACGT"*5 — contains no literal 'N'
    TagMap tags{record.Tags()};
    tags.Set(TagKey{'M', 'M'}, std::string{"N+m?,5,3;"});
    tags.Set(TagKey{'M', 'L'}, MakeUInt8Array({180, 120}));
    record.Tags(std::move(tags));

    ClipToQueryLossless(record, 4, 16);
    ASSERT_TRUE(RestoreFromLossless(record));

    const TagValue* mm{record.Tags().Get(TagKey{'M', 'M'})};
    ASSERT_NE(mm, nullptr);
    EXPECT_EQ(std::get<std::string>(*mm), "N+m?,5,3;");
    EXPECT_EQ(TagValues(record, TagKey{'M', 'L'}), (std::vector<std::uint8_t>{180, 120}));
}

TEST(LosslessClipping, RejectsUnsupportedClipTags)
{
    BamRecord record{MakeUnmappedRecord(20)};
    TagMap tags{record.Tags()};
    tags.Set(TagKey{'s', 'a'}, MakeUInt8Array({20, 1}));  // subread-pileup, not captured
    record.Tags(std::move(tags));

    EXPECT_THROW(ClipToQueryLossless(record, 4, 16), std::runtime_error);
}

TEST(LosslessClipping, NoStorageWhenNothingClipped)
{
    BamRecord record{MakeUnmappedRecord(20)};
    ClipToQueryLossless(record, 0, 20);  // full span -> no-op
    EXPECT_FALSE(record.Tags().Contains(TagKey{'l', 's'}));
    EXPECT_EQ(std::ssize(record.Sequence()), 20);
    EXPECT_FALSE(RestoreFromLossless(record));
}

}  // namespace Samoa
}  // namespace PacBio
