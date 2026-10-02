#include <pbsamoa/core/LosslessClipping.hpp>

#include <pbsamoa/core/BamRecord.hpp>
#include <pbsamoa/core/Tags.hpp>

#include <pbcopper/json/JSON.h>

#include <gtest/gtest.h>

#include <algorithm>
#include <ranges>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include <cstddef>
#include <cstdint>
#include <cstring>

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

TagArray MakeUInt32Array(const std::vector<std::uint32_t>& values)
{
    TagArray array{'I'};
    for (const std::uint32_t value : values) {
        array.AppendUInt32(value);
    }
    return array;
}

std::vector<std::uint32_t> DecodeI(const TagArray& array)
{
    std::vector<std::uint32_t> out(array.Count());
    std::memcpy(std::data(out), std::data(array.Data()), std::size(array.Data()));
    return out;
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

// The `ls` blob has lima's shape, which lima-undo depends on: both flanks list every
// modification type (an empty list when the flank holds none of its sites) under MM and ML,
// and pMM lists every type, with 0 for a type that retains no site. lima-undo rewrites MM
// sorted by type and reorders ML to match only when a flank carries an ML map, so without
// the maps a record whose MM is not sorted (as here) gets its probabilities shuffled.
// Sequence: C at 0, 2, 5, 7, 18 and A at 1, 3. "C+m?,0,1,0,0" marks the C at 0, 5, 7, 18
// (skipping the one at 2), "A+a.,0" the A at 1, and "G-h?" has no site. Clip keeps [4,16).
TEST(LosslessClipping, BaseModBlobListsEveryTypeInBothFlanksLikeLima)
{
    const std::string seq{"CACATCTCTTTTTTTTTTCT"};
    BamRecord record;
    record.Name("movie/1/ccs").Flag(0x4).Sequence(seq).Qualities(Ramp(20, 20));
    TagMap tags;
    tags.Set(TagKey{'M', 'M'}, std::string{"C+m?,0,1,0,0;A+a.,0;G-h?;"});
    tags.Set(TagKey{'M', 'L'}, MakeUInt8Array({10, 20, 30, 40, 50}));
    record.Tags(std::move(tags));

    ClipToQueryLossless(record, 4, 16);

    EXPECT_EQ(std::get<std::string>(*record.Tags().Get(TagKey{'M', 'M'})), "C+m?,0,0;A+a.;G-h?;");
    EXPECT_EQ(TagValues(record, TagKey{'M', 'L'}), (std::vector<std::uint8_t>{20, 30}));

    const JSON::Json blob = JSON::Json::from_msgpack(TagValues(record, TagKey{'l', 's'}));
    const JSON::Json expectedLeadMM = {{"C+m?", {0}}, {"A+a.", {0}}, {"G-h?", JSON::Json::array()}};
    const JSON::Json expectedLeadML = {
        {"C+m?", {10}}, {"A+a.", {50}}, {"G-h?", JSON::Json::array()}};
    const JSON::Json expectedTrailMM = {
        {"C+m?", {0}}, {"A+a.", JSON::Json::array()}, {"G-h?", JSON::Json::array()}};
    const JSON::Json expectedTrailML = {
        {"C+m?", {40}}, {"A+a.", JSON::Json::array()}, {"G-h?", JSON::Json::array()}};
    // The C at 2 precedes the first retained site; the A at 3 has no retained site to join.
    const JSON::Json expectedPrefixLost = {{"C+m?", 1}, {"A+a.", 0}, {"G-h?", 0}};
    EXPECT_EQ(blob["lead"]["MM"], expectedLeadMM);
    EXPECT_EQ(blob["lead"]["ML"], expectedLeadML);
    EXPECT_EQ(blob["trail"]["MM"], expectedTrailMM);
    EXPECT_EQ(blob["trail"]["ML"], expectedTrailML);
    EXPECT_EQ(blob["pMM"], expectedPrefixLost);

    // Restore, like lima-undo, writes the types sorted by name, with ML following.
    ASSERT_TRUE(RestoreFromLossless(record));
    EXPECT_EQ(record.Sequence(), seq);
    EXPECT_EQ(std::get<std::string>(*record.Tags().Get(TagKey{'M', 'M'})),
              "A+a.,0;C+m?,0,1,0,0;G-h?;");
    EXPECT_EQ(TagValues(record, TagKey{'M', 'L'}), (std::vector<std::uint8_t>{50, 10, 20, 30, 40}));
}

// lima notes the barcode call that led to a clip in the blob and writes the pileup fixups
// for every record, 0 without an `sa` tag. A record clipped before it is barcoded has no
// call to note.
TEST(LosslessClipping, BlobNotesBarcodeCallAndPileupFixupsLikeLima)
{
    BamRecord barcoded{MakeUnmappedRecord(20)};
    TagArray call{'S'};
    call.AppendUInt16(5);
    call.AppendUInt16(270);
    barcoded.MutableTags().Set(TagKey{'b', 'c'}, std::move(call));
    barcoded.MutableTags().Set(TagKey{'b', 'q'}, std::int64_t{87});

    ClipToQueryLossless(barcoded, 4, 16);

    const JSON::Json blob = JSON::Json::from_msgpack(TagValues(barcoded, TagKey{'l', 's'}));
    EXPECT_EQ(blob["bc"], "5--270");
    EXPECT_EQ(blob["bq"], 87);
    EXPECT_EQ(blob["psa"], 0);
    EXPECT_EQ(blob["ssa"], 0);
    EXPECT_EQ(blob["csa"], -1);
    ASSERT_TRUE(RestoreFromLossless(barcoded));

    BamRecord plain{MakeUnmappedRecord(20)};
    ClipToQueryLossless(plain, 4, 16);
    const JSON::Json plainBlob = JSON::Json::from_msgpack(TagValues(plain, TagKey{'l', 's'}));
    EXPECT_FALSE(plainBlob.contains("bc"));
    EXPECT_FALSE(plainBlob.contains("bq"));
    EXPECT_EQ(plainBlob["csa"], -1);
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
    tags.Set(TagKey{'p', 'c'}, std::string(20, 'A'));  // raw-subread pulse, not captured
    record.Tags(std::move(tags));

    EXPECT_THROW(ClipToQueryLossless(record, 4, 16), std::runtime_error);
}

// Subread pileup: sa is run-length encoded, so a clip boundary can fall inside a run.
// The straddling runs must shrink on clip and grow back on undo (lima psa/ssa), and the
// flank-only runs must move out to `ls` and back.
TEST(LosslessClipping, RoundTripRestoresSubreadPileupTags)
{
    BamRecord record{MakeUnmappedRecord(20)};
    TagMap tags{record.Tags()};
    tags.Set(TagKey{'s', 'a'}, MakeUInt32Array({2, 5, 6, 7, 8, 9, 4, 3}));  // runs: 2+6+8+4
    tags.Set(TagKey{'s', 'm'}, MakeUInt8Array(Ramp(200, 20)));
    tags.Set(TagKey{'s', 'x'}, MakeUInt8Array(Ramp(60, 20)));
    record.Tags(std::move(tags));
    const BamRecord original{record};

    ClipToQueryLossless(record, 4, 16);  // splits run 2 (6 -> 4) and run 3 (8 -> 8 kept, 4 lost)
    EXPECT_EQ(DecodeI(ArrayTag(record, TagKey{'s', 'a'})),
              (std::vector<std::uint32_t>{4, 7, 8, 9}));
    EXPECT_EQ(TagValues(record, TagKey{'s', 'm'}), Ramp(204, 12));
    EXPECT_EQ(TagValues(record, TagKey{'s', 'x'}), Ramp(64, 12));

    ASSERT_TRUE(RestoreFromLossless(record));
    EXPECT_EQ(DecodeI(ArrayTag(record, TagKey{'s', 'a'})),
              DecodeI(ArrayTag(original, TagKey{'s', 'a'})));
    EXPECT_EQ(TagValues(record, TagKey{'s', 'm'}), TagValues(original, TagKey{'s', 'm'}));
    EXPECT_EQ(TagValues(record, TagKey{'s', 'x'}), TagValues(original, TagKey{'s', 'x'}));
    EXPECT_EQ(record.Sequence(), original.Sequence());
}

TEST(LosslessClipping, NestedClipsRestoreSubreadPileupTags)
{
    BamRecord record{MakeUnmappedRecord(20)};
    TagMap tags{record.Tags()};
    tags.Set(TagKey{'s', 'a'}, MakeUInt32Array({2, 5, 6, 7, 8, 9, 4, 3}));
    tags.Set(TagKey{'s', 'm'}, MakeUInt8Array(Ramp(200, 20)));
    tags.Set(TagKey{'s', 'x'}, MakeUInt8Array(Ramp(60, 20)));
    record.Tags(std::move(tags));
    const BamRecord original{record};

    ClipToQueryLossless(record, 4, 16);
    ClipToQueryLossless(record, 1, 11);  // second boundary inside the same runs again
    EXPECT_EQ(DecodeI(ArrayTag(record, TagKey{'s', 'a'})),
              (std::vector<std::uint32_t>{3, 7, 7, 9}));

    ASSERT_TRUE(RestoreFromLossless(record));
    EXPECT_EQ(DecodeI(ArrayTag(record, TagKey{'s', 'a'})),
              DecodeI(ArrayTag(original, TagKey{'s', 'a'})));
    EXPECT_EQ(TagValues(record, TagKey{'s', 'm'}), TagValues(original, TagKey{'s', 'm'}));
    EXPECT_EQ(TagValues(record, TagKey{'s', 'x'}), TagValues(original, TagKey{'s', 'x'}));
}

// Undo also reverts what a demultiplexer stamps around the clip: barcode bits in cx,
// the barcoded read-group suffix, and (CCS only) the qs/qe query window.
TEST(LosslessClipping, RestoreStripsDemultiplexingSideEffects)
{
    BamRecord record{MakeUnmappedRecord(20)};
    ClipToQueryLossless(record, 4, 16);
    TagMap tags{record.Tags()};
    tags.Set(TagKey{'c', 'x'}, std::int64_t{4 | 8 | 16});  // BARCODE_BEFORE|BARCODE_AFTER|FORWARD
    tags.Set(TagKey{'q', 's'}, std::int64_t{4});
    tags.Set(TagKey{'q', 'e'}, std::int64_t{16});
    tags.Set(RG_TAG, std::string{"abcd1234/1--1"});
    tags.Set(TagKey{'b', 'c'}, MakeUInt8Array({1, 1}));
    record.Tags(std::move(tags));

    ASSERT_TRUE(RestoreFromLossless(record));
    EXPECT_EQ(std::get<std::int64_t>(*record.Tags().Get(TagKey{'c', 'x'})), 16);
    EXPECT_FALSE(record.Tags().Contains(TagKey{'q', 's'}));
    EXPECT_FALSE(record.Tags().Contains(TagKey{'q', 'e'}));
    EXPECT_FALSE(record.Tags().Contains(TagKey{'b', 'c'}));
    EXPECT_EQ(record.ReadGroupId(), "abcd1234");
}

TEST(LosslessClipping, RestoreDropsEmptyCxOnCcsOnly)
{
    BamRecord ccs{MakeUnmappedRecord(20)};
    ClipToQueryLossless(ccs, 4, 16);
    ccs.MutableTags().Set(TagKey{'c', 'x'}, std::int64_t{4 | 8});
    ASSERT_TRUE(RestoreFromLossless(ccs));
    EXPECT_FALSE(ccs.Tags().Contains(TagKey{'c', 'x'}));

    BamRecord subread{MakeUnmappedRecord(20)};
    subread.Name("movie/1/0_20");
    ClipToQueryLossless(subread, 4, 16);
    subread.MutableTags().Set(TagKey{'c', 'x'}, std::int64_t{4 | 8});
    ASSERT_TRUE(RestoreFromLossless(subread));
    EXPECT_EQ(std::get<std::int64_t>(*subread.Tags().Get(TagKey{'c', 'x'})), 0);
}

// dq/iq/mq/sq/dt/st are per-base tags the PacBio default clipper trims (substring
// strategy) but lossless storage neither captures nor restores. Clipping a record
// carrying them would silently shorten them with no way to undo, so they must be
// rejected up front like the pulse/pileup tags.
TEST(LosslessClipping, RejectsDefaultClippedPerBaseTags)
{
    for (const TagKey key : {TagKey{'d', 'q'}, TagKey{'i', 'q'}, TagKey{'m', 'q'}, TagKey{'s', 'q'},
                             TagKey{'d', 't'}, TagKey{'s', 't'}}) {
        const std::string label{key.First(), key.Second()};
        BamRecord record{MakeUnmappedRecord(20)};
        TagMap tags{record.Tags()};
        tags.Set(key, std::string(20, '!'));  // per-base length matches the sequence
        record.Tags(std::move(tags));
        EXPECT_THROW(ClipToQueryLossless(record, 4, 16), std::runtime_error) << label;
    }
}

TEST(LosslessClipping, NoStorageWhenNothingClipped)
{
    BamRecord record{MakeUnmappedRecord(20)};
    ClipToQueryLossless(record, 0, 20);  // full span -> no-op
    EXPECT_FALSE(record.Tags().Contains(TagKey{'l', 's'}));
    EXPECT_EQ(std::ssize(record.Sequence()), 20);
    EXPECT_FALSE(RestoreFromLossless(record));
}

// An ML array that does not cover every MM site cannot be sliced losslessly (undo would
// lose probabilities), so the clip must reject rather than silently drop data — the same
// contract that rejects uncaptured per-base tags.
TEST(LosslessClipping, RejectsMlShorterThanMmSites)
{
    BamRecord record{MakeUnmappedRecord(20)};
    TagMap tags{record.Tags()};
    tags.Set(TagKey{'M', 'M'}, std::string{"C+m?,0,1,1;"});  // three modification sites
    tags.Set(TagKey{'M', 'L'}, MakeUInt8Array({200, 150}));  // only two probabilities
    record.Tags(std::move(tags));

    EXPECT_THROW(ClipToQueryLossless(record, 4, 16), std::runtime_error);
}

// MN:i records the SEQ length MM/ML were produced against; the clip rewrites it to the
// clipped length and restore must return it to the full length (htslib rejects MN != l_qseq).
TEST(LosslessClipping, RoundTripRestoresMnTag)
{
    BamRecord record{MakeUnmappedRecord(20)};
    TagMap tags{record.Tags()};
    tags.Set(TagKey{'M', 'M'}, std::string{"C+m?,0,1,1;"});
    tags.Set(TagKey{'M', 'L'}, MakeUInt8Array({200, 150, 100}));
    tags.Set(TagKey{'M', 'N'}, std::int64_t{20});
    record.Tags(std::move(tags));

    ClipToQueryLossless(record, 4, 16);  // 20 -> 12
    const TagValue* clippedMn{record.Tags().Get(TagKey{'M', 'N'})};
    ASSERT_NE(clippedMn, nullptr);
    EXPECT_EQ(std::get<std::int64_t>(*clippedMn), 12);

    ASSERT_TRUE(RestoreFromLossless(record));
    const TagValue* restoredMn{record.Tags().Get(TagKey{'M', 'N'})};
    ASSERT_NE(restoredMn, nullptr);
    EXPECT_EQ(std::get<std::int64_t>(*restoredMn), 20);
}

}  // namespace Samoa
}  // namespace PacBio
