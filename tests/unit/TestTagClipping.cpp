#include <pbsamoa/core/TagClipping.hpp>
#include <pbsamoa/core/Tags.hpp>

#include <gtest/gtest.h>

#include <string>
#include <vector>

#include <cstddef>
#include <cstdint>
#include <cstring>

namespace PacBio {
namespace Samoa {

namespace {
TagArray MakeUInt8Array(std::initializer_list<std::uint8_t> values)
{
    TagArray arr{'C'};
    for (const auto v : values) {
        arr.AppendUInt8(v);
    }
    return arr;
}

TagArray MakeUInt16Array(std::initializer_list<std::uint16_t> values)
{
    TagArray arr{'S'};
    for (const auto v : values) {
        arr.AppendUInt16(v);
    }
    return arr;
}
}  // namespace

TEST(TagClipping, SubstringClip_String_TrimBothEnds)
{
    const SubstringClipStrategy strategy;
    TagValue val{std::string{"ABCDEFGHIJ"}};
    ASSERT_TRUE(strategy.Clip(val, 2, 6, 10));
    EXPECT_EQ(std::get<std::string>(val), "CDEFGH");
}

TEST(TagClipping, SubstringClip_String_NoTrim)
{
    const SubstringClipStrategy strategy;
    TagValue val{std::string{"ABCDE"}};
    ASSERT_TRUE(strategy.Clip(val, 0, 5, 5));
    EXPECT_EQ(std::get<std::string>(val), "ABCDE");
}

TEST(TagClipping, SubstringClip_String_Empty)
{
    const SubstringClipStrategy strategy;
    TagValue val{std::string{"ABC"}};
    ASSERT_TRUE(strategy.Clip(val, 1, 0, 3));
    EXPECT_EQ(std::get<std::string>(val), "");
}

TEST(TagClipping, SubstringClip_UInt8Array_TrimFront)
{
    const SubstringClipStrategy strategy;
    TagValue val{MakeUInt8Array({10, 20, 30, 40, 50})};
    ASSERT_TRUE(strategy.Clip(val, 2, 3, 5));
    const auto& arr{std::get<TagArray>(val)};
    EXPECT_EQ(arr.Count(), 3U);
    EXPECT_EQ(static_cast<std::uint8_t>(arr.Data()[0]), 30);
    EXPECT_EQ(static_cast<std::uint8_t>(arr.Data()[1]), 40);
    EXPECT_EQ(static_cast<std::uint8_t>(arr.Data()[2]), 50);
}

TEST(TagClipping, SubstringClip_UInt16Array)
{
    const SubstringClipStrategy strategy;
    TagValue val{MakeUInt16Array({100, 200, 300, 400, 500})};
    ASSERT_TRUE(strategy.Clip(val, 1, 3, 5));
    const auto& arr{std::get<TagArray>(val)};
    EXPECT_EQ(arr.Count(), 3U);
    std::uint16_t v0;
    std::memcpy(&v0, std::data(arr.Data()), sizeof(v0));
    EXPECT_EQ(v0, 200);
}

TEST(TagClipping, SubstringClip_UnsupportedType_ReturnsFalse)
{
    const SubstringClipStrategy strategy;
    TagValue val{std::int64_t{42}};
    EXPECT_FALSE(strategy.Clip(val, 0, 1, 1));
}

TEST(TagClipping, ReverseSubstringClip_String)
{
    const ReverseSubstringClipStrategy strategy;
    // clipOffset=1, clipLength=3, seqLen=10
    // reverseOffset = 10 - (1+3) = 6
    TagValue val{std::string{"0123456789"}};
    ASSERT_TRUE(strategy.Clip(val, 1, 3, 10));
    EXPECT_EQ(std::get<std::string>(val), "678");
}

TEST(TagClipping, ReverseSubstringClip_UInt8Array)
{
    const ReverseSubstringClipStrategy strategy;
    // clipOffset=0, clipLength=2, seqLen=5
    // reverseOffset = 5 - (0+2) = 3
    TagValue val{MakeUInt8Array({10, 20, 30, 40, 50})};
    ASSERT_TRUE(strategy.Clip(val, 0, 2, 5));
    const auto& arr{std::get<TagArray>(val)};
    EXPECT_EQ(arr.Count(), 2U);
    EXPECT_EQ(static_cast<std::uint8_t>(arr.Data()[0]), 40);
    EXPECT_EQ(static_cast<std::uint8_t>(arr.Data()[1]), 50);
}

TEST(TagClipping, TagClipper_ClipsRegisteredTags)
{
    TagClipper clipper;
    static const SubstringClipStrategy strategy;
    clipper.Register({TagKey{'i', 'p'}, TagKey{'p', 'w'}}, strategy);

    TagMap tags;
    tags.Set(TagKey{'i', 'p'}, MakeUInt8Array({10, 20, 30, 40, 50}));
    tags.Set(TagKey{'p', 'w'}, MakeUInt8Array({1, 2, 3, 4, 5}));
    tags.Set(TagKey{'N', 'M'}, std::int64_t{42});  // unregistered

    clipper.ClipTags(tags, 1, 3, 5);

    const auto* ip{tags.Get(TagKey{'i', 'p'})};
    ASSERT_NE(ip, nullptr);
    EXPECT_EQ(std::get<TagArray>(*ip).Count(), 3U);

    const auto* pw{tags.Get(TagKey{'p', 'w'})};
    ASSERT_NE(pw, nullptr);
    EXPECT_EQ(std::get<TagArray>(*pw).Count(), 3U);

    // NM untouched
    const auto* nm{tags.Get(TagKey{'N', 'M'})};
    ASSERT_NE(nm, nullptr);
    EXPECT_EQ(std::get<std::int64_t>(*nm), 42);
}

TEST(TagClipping, TagClipper_MixedStrategies)
{
    TagClipper clipper;
    static const SubstringClipStrategy substringStrategy;
    clipper.Register({TagKey{'i', 'p'}}, substringStrategy);
    static const ReverseSubstringClipStrategy reverseSubstringStrategy;
    clipper.Register({TagKey{'r', 'i'}}, reverseSubstringStrategy);

    TagMap tags;
    tags.Set(TagKey{'i', 'p'}, MakeUInt8Array({10, 20, 30, 40, 50}));
    tags.Set(TagKey{'r', 'i'}, MakeUInt8Array({50, 40, 30, 20, 10}));

    // clipOffset=1, clipLength=3, seqLen=5
    clipper.ClipTags(tags, 1, 3, 5);

    // ip: substring [1,4) => {20, 30, 40}
    const auto* ip{tags.Get(TagKey{'i', 'p'})};
    ASSERT_NE(ip, nullptr);
    EXPECT_EQ(std::get<TagArray>(*ip).Count(), 3U);
    EXPECT_EQ(static_cast<std::uint8_t>(std::get<TagArray>(*ip).Data()[0]), 20);

    // ri: reverse offset = 5 - (1+3) = 1, substring [1,4) => {40, 30, 20}
    const auto* ri{tags.Get(TagKey{'r', 'i'})};
    ASSERT_NE(ri, nullptr);
    EXPECT_EQ(std::get<TagArray>(*ri).Count(), 3U);
    EXPECT_EQ(static_cast<std::uint8_t>(std::get<TagArray>(*ri).Data()[0]), 40);
}

TEST(TagClipping, TagClipper_RemovesTagOnClipFailure)
{
    TagClipper clipper;
    static const SubstringClipStrategy baStrategy;
    clipper.Register({TagKey{'b', 'a'}}, baStrategy);

    TagMap tags;
    // int64 is unsupported by SubstringClipStrategy => Clip returns false => tag removed
    tags.Set(TagKey{'b', 'a'}, std::int64_t{99});

    clipper.ClipTags(tags, 0, 1, 1);

    EXPECT_FALSE(tags.Contains(TagKey{'b', 'a'}));
}

TEST(TagClipping, PacBioDefault_RegistersExpectedTags)
{
    const TagClipper clipper{TagClipper::PacBioDefault()};

    // Build a TagMap with all expected PacBio tags as uint8 arrays of length 5
    TagMap tags;
    const std::vector<TagKey> substringKeys{
        TagKey{'d', 'q'}, TagKey{'i', 'q'}, TagKey{'m', 'q'}, TagKey{'s', 'q'}, TagKey{'d', 't'},
        TagKey{'s', 't'}, TagKey{'i', 'p'}, TagKey{'p', 'w'}, TagKey{'f', 'i'}, TagKey{'f', 'p'},
    };
    const std::vector<TagKey> reverseKeys{
        TagKey{'r', 'i'},
        TagKey{'r', 'p'},
    };

    for (const auto& key : substringKeys) {
        tags.Set(key, MakeUInt8Array({10, 20, 30, 40, 50}));
    }
    for (const auto& key : reverseKeys) {
        tags.Set(key, MakeUInt8Array({10, 20, 30, 40, 50}));
    }

    // clipOffset=1, clipLength=3, seqLen=5
    clipper.ClipTags(tags, 1, 3, 5);

    // Substring tags: [1,4) => {20, 30, 40}
    for (const auto& key : substringKeys) {
        const auto* val{tags.Get(key)};
        ASSERT_NE(val, nullptr) << "missing: " << key.First() << key.Second();
        EXPECT_EQ(std::get<TagArray>(*val).Count(), 3U);
        EXPECT_EQ(static_cast<std::uint8_t>(std::get<TagArray>(*val).Data()[0]), 20);
    }

    // Reverse tags: reverseOffset = 5 - (1+3) = 1, [1,4) => {20, 30, 40}
    for (const auto& key : reverseKeys) {
        const auto* val{tags.Get(key)};
        ASSERT_NE(val, nullptr) << "missing: " << key.First() << key.Second();
        EXPECT_EQ(std::get<TagArray>(*val).Count(), 3U);
        EXPECT_EQ(static_cast<std::uint8_t>(std::get<TagArray>(*val).Data()[0]), 20);
    }
}

// --- PulseClipStrategy tests ---

TEST(TagClipping, PulseClip_BasicMapping)
{
    // pc = "AaCgTt" (6 pulses, basecalled at 0, 2, 4 = bases A, C, T)
    // 3 bases total. Clip bases [1, 2) => keep just base index 1 (C at pulse 2)
    // Pulse range: startPulse = FindNthBase(1) = 2, endPulse = FindNthBase(1) = 2
    // Pulse clip: [2, 3) = single pulse "C"

    TagMap tags;
    tags.Set(TagKey{'p', 'c'}, std::string{"AaCgTt"});
    // pt is a companion string pulse tag
    tags.Set(TagKey{'p', 't'}, std::string{"123456"});

    TagClipper clipper;
    static const PulseClipStrategy pulseStrategy1;
    clipper.Register({TagKey{'p', 'c'}, TagKey{'p', 't'}}, pulseStrategy1);

    clipper.ClipTags(tags, 1, 1, 3);

    const auto* pc{tags.Get(TagKey{'p', 'c'})};
    ASSERT_NE(pc, nullptr);
    EXPECT_EQ(std::get<std::string>(*pc), "C");

    const auto* pt{tags.Get(TagKey{'p', 't'})};
    ASSERT_NE(pt, nullptr);
    EXPECT_EQ(std::get<std::string>(*pt), "3");
}

TEST(TagClipping, PulseClip_KeepMultipleBases)
{
    // pc = "AaCgTt" (basecalled: 0=A, 2=C, 4=T)
    // Clip bases [0, 2) => keep base 0 (A, pulse 0) and base 1 (C, pulse 2)
    // startPulse=0, endPulse=2 => pulse range [0,3) = "AaC"
    TagMap tags;
    tags.Set(TagKey{'p', 'c'}, std::string{"AaCgTt"});
    tags.Set(TagKey{'p', 't'}, std::string{"XYZWUV"});

    TagClipper clipper;
    static const PulseClipStrategy pulseStrategy2;
    clipper.Register({TagKey{'p', 'c'}, TagKey{'p', 't'}}, pulseStrategy2);

    clipper.ClipTags(tags, 0, 2, 3);

    EXPECT_EQ(std::get<std::string>(*tags.Get(TagKey{'p', 'c'})), "AaC");
    EXPECT_EQ(std::get<std::string>(*tags.Get(TagKey{'p', 't'})), "XYZ");
}

TEST(TagClipping, PulseClip_AllBases)
{
    // pc = "AaCgTt", clip all 3 bases [0,3)
    // startPulse=0, endPulse=4 => range [0, 5) = "AaCgT"
    // The trailing "t" after the last basecalled T is NOT included.
    TagMap tags;
    tags.Set(TagKey{'p', 'c'}, std::string{"AaCgTt"});

    TagClipper clipper;
    static const PulseClipStrategy pulseStrategy3;
    clipper.Register({TagKey{'p', 'c'}}, pulseStrategy3);

    clipper.ClipTags(tags, 0, 3, 3);

    EXPECT_EQ(std::get<std::string>(*tags.Get(TagKey{'p', 'c'})), "AaCgT");
}

TEST(TagClipping, PulseClip_UInt16Array)
{
    // pc = "AaCT" (4 pulses, basecalled: 0=A, 2=C, 3=T)
    // pa is a uint16 pulse tag
    // Clip bases [1, 3) => base 1 (C, pulse 2) and base 2 (T, pulse 3)
    // startPulse=2, endPulse=3 => pulse range [2,4)

    TagMap tags;
    tags.Set(TagKey{'p', 'c'}, std::string{"AaCT"});
    TagArray pa{'S'};
    pa.AppendUInt16(100);
    pa.AppendUInt16(200);
    pa.AppendUInt16(300);
    pa.AppendUInt16(400);
    tags.Set(TagKey{'p', 'a'}, std::move(pa));

    TagClipper clipper;
    static const PulseClipStrategy pulseStrategy4;
    clipper.Register({TagKey{'p', 'c'}, TagKey{'p', 'a'}}, pulseStrategy4);

    clipper.ClipTags(tags, 1, 2, 3);

    const auto& paResult{std::get<TagArray>(*tags.Get(TagKey{'p', 'a'}))};
    EXPECT_EQ(paResult.Count(), 2U);
    std::uint16_t v0;
    std::memcpy(&v0, std::data(paResult.Data()), sizeof(v0));
    EXPECT_EQ(v0, 300);
}

// --- BasemodClipStrategy tests ---

TEST(TagClipping, BasemodClip_SimpleClip)
{
    // Sequence: "AACAACAA" (8 bases, C at positions 2 and 5)
    // MM: "C+m,0,0;" means:
    //   Walk C's in sequence: C[0] at pos 2, C[1] at pos 5
    //   Skip 0 -> C[0] (pos 2) is modified
    //   Skip 0 -> C[1] (pos 5) is modified
    // ML: [200, 100] (quality per mod site)
    //
    // Clip to [3, 8) => keep "AACAA" (positions 3..7)
    // In the clipped window: only C at position 5 exists
    // Only 1 C is before clip window (pos 2). basesBeforeClip=1, basesInClip=1
    // prefixSum = {1, 2}. lower_bound(2) -> index 1. upper_bound(2) -> index 2.
    // frontRemoved=1, retained=1. newSkips[0] = prefixSum[1] - 1 - 1 = 0
    // New MM: "C+m,0;" New ML: [100]

    TagMap tags;
    tags.Set(TagKey{'M', 'M'}, std::string{"C+m,0,0;"});
    tags.Set(TagKey{'M', 'L'}, MakeUInt8Array({200, 100}));

    TagClipper clipper;
    static const BasemodClipStrategy basemodStrategy;
    clipper.Register({TagKey{'M', 'M'}, TagKey{'M', 'L'}}, basemodStrategy);

    clipper.ClipTags(tags, 3, 5, 8, "AACAACAA");

    const auto* mm{tags.Get(TagKey{'M', 'M'})};
    ASSERT_NE(mm, nullptr);
    EXPECT_EQ(std::get<std::string>(*mm), "C+m,0;");

    const auto* ml{tags.Get(TagKey{'M', 'L'})};
    ASSERT_NE(ml, nullptr);
    EXPECT_EQ(std::get<TagArray>(*ml).Count(), 1U);
    EXPECT_EQ(static_cast<std::uint8_t>(std::get<TagArray>(*ml).Data()[0]), 100);
}

TEST(TagClipping, BasemodClip_NoModsInWindow)
{
    // Sequence: "AACAATAA" (C at position 2 only; no C in [3,8))
    // MM: "C+m,0;" => C[0] at pos 2 is modified
    // ML: [200]
    // Clip to [3, 5) => "AA" — no C's in window
    // basesBeforeClip=1, basesInClip=0
    // prefixSum = {1}. lower_bound(2) -> end. upper_bound(1) -> end.
    // frontRemoved=1, retained=0
    // New MM: "C+m;" (type header, no mods)
    // New ML: empty

    TagMap tags;
    tags.Set(TagKey{'M', 'M'}, std::string{"C+m,0;"});
    tags.Set(TagKey{'M', 'L'}, MakeUInt8Array({200}));

    TagClipper clipper;
    static const BasemodClipStrategy basemodStrategy;
    clipper.Register({TagKey{'M', 'M'}, TagKey{'M', 'L'}}, basemodStrategy);

    clipper.ClipTags(tags, 3, 5, 8, "AACAATAA");

    const auto* mm{tags.Get(TagKey{'M', 'M'})};
    ASSERT_NE(mm, nullptr);
    EXPECT_EQ(std::get<std::string>(*mm), "C+m;");

    const auto* ml{tags.Get(TagKey{'M', 'L'})};
    ASSERT_NE(ml, nullptr);
    EXPECT_EQ(std::get<TagArray>(*ml).Count(), 0U);
}

TEST(TagClipping, BasemodClip_AllModsRetained)
{
    // Sequence: "CCCC" (all C's)
    // MM: "C+m,0,1;" => C[0] modified, skip 1, C[2] modified
    // ML: [200, 100]
    // Clip to [0, 4) => full sequence, no change
    // basesBeforeClip=0, basesInClip=4
    // prefixSum = {1, 3}. lower_bound(1) -> index 0. upper_bound(4) -> end.
    // frontRemoved=0, retained=2
    // newSkips = {0, 1}, skip[0] adjusted: prefixSum[0] - 0 - 1 = 0 (unchanged)
    // New MM: "C+m,0,1;" New ML: [200, 100]

    TagMap tags;
    tags.Set(TagKey{'M', 'M'}, std::string{"C+m,0,1;"});
    tags.Set(TagKey{'M', 'L'}, MakeUInt8Array({200, 100}));

    TagClipper clipper;
    static const BasemodClipStrategy basemodStrategy;
    clipper.Register({TagKey{'M', 'M'}, TagKey{'M', 'L'}}, basemodStrategy);

    clipper.ClipTags(tags, 0, 4, 4, "CCCC");

    EXPECT_EQ(std::get<std::string>(*tags.Get(TagKey{'M', 'M'})), "C+m,0,1;");
    const auto& ml{std::get<TagArray>(*tags.Get(TagKey{'M', 'L'}))};
    EXPECT_EQ(ml.Count(), 2U);
    EXPECT_EQ(static_cast<std::uint8_t>(ml.Data()[0]), 200);
    EXPECT_EQ(static_cast<std::uint8_t>(ml.Data()[1]), 100);
}

TEST(TagClipping, BasemodClip_WithQuestionMark)
{
    // Sequence: "AACCGATCC" (C at positions 2,3,7,8)
    // MM: "C+m?,1,0;" => skip 1 C, mod C[1](pos 3); skip 0, mod C[2](pos 7)
    // ML: [220, 180]
    // Clip to [3, 6) => keep "CGA" (positions 3..5), C at position 3 only
    // basesBeforeClip=2 (C at pos 2 and 3... wait, pos 3 is within [3,6))
    // Actually: bases before clip (positions 0..2): C at pos 2 -> 1 C
    // Bases in clip (positions 3..5): C at pos 3 -> 1 C
    // prefixSum = {2, 3}. lower_bound(2) -> index 0. upper_bound(2) -> index 1.
    // frontRemoved=0, retained=1
    // newSkips[0] = prefixSum[0] - 1 - 1 = 0
    // New MM: "C+m?,0;" New ML: [220]

    TagMap tags;
    tags.Set(TagKey{'M', 'M'}, std::string{"C+m?,1,0;"});
    tags.Set(TagKey{'M', 'L'}, MakeUInt8Array({220, 180}));

    TagClipper clipper;
    static const BasemodClipStrategy basemodStrategy;
    clipper.Register({TagKey{'M', 'M'}, TagKey{'M', 'L'}}, basemodStrategy);

    clipper.ClipTags(tags, 3, 3, 9, "AACCGATCC");

    EXPECT_EQ(std::get<std::string>(*tags.Get(TagKey{'M', 'M'})), "C+m?,0;");
    const auto& ml{std::get<TagArray>(*tags.Get(TagKey{'M', 'L'}))};
    EXPECT_EQ(ml.Count(), 1U);
    EXPECT_EQ(static_cast<std::uint8_t>(ml.Data()[0]), 220);
}

TEST(TagClipping, BasemodClip_MultipleModTypes)
{
    // Sequence: "ACGTACGT" (8 bases)
    // C at positions 1, 5; A at positions 0, 4
    // MM: "C+m,0,0;A+a,0,0;" => both C's and both A's modified
    // ML: [200, 100, 150, 250] (C+m QVs first, then A+a QVs)
    //
    // Clip to [4, 4) => "ACGT" (positions 4..7)
    // C+m: frontRemoved=1 (pos 1), retained=1 (pos 5) -> keep ML[1]=100
    // A+a: frontRemoved=1 (pos 0), retained=1 (pos 4) -> keep ML[3]=250
    // New MM: "C+m,0;A+a,0;" New ML: [100, 250]
    //
    // This exercises per-type ML slicing (not a flat window).

    TagMap tags;
    tags.Set(TagKey{'M', 'M'}, std::string{"C+m,0,0;A+a,0,0;"});
    tags.Set(TagKey{'M', 'L'}, MakeUInt8Array({200, 100, 150, 250}));

    TagClipper clipper;
    static const BasemodClipStrategy basemodStrategy;
    clipper.Register({TagKey{'M', 'M'}, TagKey{'M', 'L'}}, basemodStrategy);

    clipper.ClipTags(tags, 4, 4, 8, "ACGTACGT");

    EXPECT_EQ(std::get<std::string>(*tags.Get(TagKey{'M', 'M'})), "C+m,0;A+a,0;");
    const auto& ml{std::get<TagArray>(*tags.Get(TagKey{'M', 'L'}))};
    EXPECT_EQ(ml.Count(), 2U);
    EXPECT_EQ(static_cast<std::uint8_t>(ml.Data()[0]), 100);
    EXPECT_EQ(static_cast<std::uint8_t>(ml.Data()[1]), 250);
}

TEST(TagClipping, BasemodClip_EmptySequence_RemovesTags)
{
    // When no sequence is provided, basemods tags should be removed
    TagMap tags;
    tags.Set(TagKey{'M', 'M'}, std::string{"C+m,0;"});
    tags.Set(TagKey{'M', 'L'}, MakeUInt8Array({200}));

    TagClipper clipper;
    static const BasemodClipStrategy basemodStrategy;
    clipper.Register({TagKey{'M', 'M'}, TagKey{'M', 'L'}}, basemodStrategy);

    // No sequence provided (empty string_view)
    clipper.ClipTags(tags, 0, 4, 4);

    EXPECT_FALSE(tags.Contains(TagKey{'M', 'M'}));
    EXPECT_FALSE(tags.Contains(TagKey{'M', 'L'}));
}

TEST(TagClipping, BasemodClip_TrimFrontAndBack)
{
    // Sequence: "CCCCC" (5 C's at positions 0..4)
    // MM: "C+m,0,0,0,0,0;" => all 5 C's modified
    // ML: [10, 20, 30, 40, 50]
    // Clip to [1, 3) => keep "CCC" (positions 1..3)
    // basesBeforeClip=1, basesInClip=3
    // prefixSum = {1, 2, 3, 4, 5}
    // lower_bound(2) -> idx 1. upper_bound(4) -> idx 4.
    // frontRemoved=1, retained=3
    // newSkips = {0, 0, 0} -> skip[0] adjusted: prefixSum[1] - 1 - 1 = 0 (no change)
    // New MM: "C+m,0,0,0;" New ML: [20, 30, 40]

    TagMap tags;
    tags.Set(TagKey{'M', 'M'}, std::string{"C+m,0,0,0,0,0;"});
    tags.Set(TagKey{'M', 'L'}, MakeUInt8Array({10, 20, 30, 40, 50}));

    TagClipper clipper;
    static const BasemodClipStrategy basemodStrategy;
    clipper.Register({TagKey{'M', 'M'}, TagKey{'M', 'L'}}, basemodStrategy);

    clipper.ClipTags(tags, 1, 3, 5, "CCCCC");

    EXPECT_EQ(std::get<std::string>(*tags.Get(TagKey{'M', 'M'})), "C+m,0,0,0;");
    const auto& ml{std::get<TagArray>(*tags.Get(TagKey{'M', 'L'}))};
    EXPECT_EQ(ml.Count(), 3U);
    EXPECT_EQ(static_cast<std::uint8_t>(ml.Data()[0]), 20);
    EXPECT_EQ(static_cast<std::uint8_t>(ml.Data()[1]), 30);
    EXPECT_EQ(static_cast<std::uint8_t>(ml.Data()[2]), 40);
}

TEST(TagClipping, BasemodClip_SkipAdjustment)
{
    // Sequence: "ATCATCATC" (9 bases, C at positions 2, 5, 8)
    // MM: "C+m,2;" => skip 2 C's, modify C[2] (pos 8)
    // ML: [255]
    // Clip to [4, 5) => keep "TCATC" (positions 4..8)
    // basesBeforeClip(C)=1 (C at pos 2). basesInClip(C)=2 (C at pos 5, 8)
    // prefixSum = {3}. lower_bound(2) -> idx 0. upper_bound(3) -> idx 1.
    // frontRemoved=0, retained=1
    // newSkips[0] = prefixSum[0] - 1 - 1 = 1
    // New MM: "C+m,1;" (skip 1 C in the clip window, then modify the 2nd)
    // New ML: [255]

    TagMap tags;
    tags.Set(TagKey{'M', 'M'}, std::string{"C+m,2;"});
    tags.Set(TagKey{'M', 'L'}, MakeUInt8Array({255}));

    TagClipper clipper;
    static const BasemodClipStrategy basemodStrategy;
    clipper.Register({TagKey{'M', 'M'}, TagKey{'M', 'L'}}, basemodStrategy);

    clipper.ClipTags(tags, 4, 5, 9, "ATCATCATC");

    EXPECT_EQ(std::get<std::string>(*tags.Get(TagKey{'M', 'M'})), "C+m,1;");
    const auto& ml{std::get<TagArray>(*tags.Get(TagKey{'M', 'L'}))};
    EXPECT_EQ(ml.Count(), 1U);
    EXPECT_EQ(static_cast<std::uint8_t>(ml.Data()[0]), 255);
}

// --- PileupClipStrategy tests ---

TEST(TagClipping, PileupClip_Sa_SimpleClip)
{
    // sa = [4, 10, 3, 20, 3, 15] meaning:
    //   positions 0-3: coverage 10 (run of 4)
    //   positions 4-6: coverage 20 (run of 3)
    //   positions 7-9: coverage 15 (run of 3)
    // Total length = 10
    // Clip to [2, 8) = positions 2-7, clipLength=6
    // Result:
    //   positions 2-3 from run 1: length 2, cov 10
    //   positions 4-6 from run 2: length 3, cov 20
    //   position 7 from run 3: length 1, cov 15
    // sa = [2, 10, 3, 20, 1, 15]

    const PileupClipStrategy strategy;
    TagValue val{MakeUInt8Array({4, 10, 3, 20, 3, 15})};

    ASSERT_TRUE(strategy.Clip(val, 2, 6, 10));
    const auto& result{std::get<TagArray>(val)};
    EXPECT_EQ(result.Count(), 6U);  // 3 pairs = 6 bytes
    const auto data{result.Data()};
    EXPECT_EQ(static_cast<std::uint8_t>(data[0]), 2);  // trimmed run 1
    EXPECT_EQ(static_cast<std::uint8_t>(data[1]), 10);
    EXPECT_EQ(static_cast<std::uint8_t>(data[2]), 3);  // full run 2
    EXPECT_EQ(static_cast<std::uint8_t>(data[3]), 20);
    EXPECT_EQ(static_cast<std::uint8_t>(data[4]), 1);  // trimmed run 3
    EXPECT_EQ(static_cast<std::uint8_t>(data[5]), 15);
}

TEST(TagClipping, PileupClip_Sa_WholeRuns)
{
    // sa = [3, 10, 3, 20, 4, 15] total = 10
    // Clip to [3, 6) = exactly run 2
    // Result: [3, 20]
    const PileupClipStrategy strategy;
    TagValue val{MakeUInt8Array({3, 10, 3, 20, 4, 15})};

    ASSERT_TRUE(strategy.Clip(val, 3, 3, 10));
    const auto& result{std::get<TagArray>(val)};
    EXPECT_EQ(result.Count(), 2U);
    const auto data{result.Data()};
    EXPECT_EQ(static_cast<std::uint8_t>(data[0]), 3);
    EXPECT_EQ(static_cast<std::uint8_t>(data[1]), 20);
}

TEST(TagClipping, PileupClip_Sa_SingleRun)
{
    // sa = [10, 5] (single run)
    // Clip to [2, 7) => [5, 5]
    const PileupClipStrategy strategy;
    TagValue val{MakeUInt8Array({10, 5})};

    ASSERT_TRUE(strategy.Clip(val, 2, 5, 10));
    const auto& result{std::get<TagArray>(val)};
    EXPECT_EQ(result.Count(), 2U);
    const auto data{result.Data()};
    EXPECT_EQ(static_cast<std::uint8_t>(data[0]), 5);
    EXPECT_EQ(static_cast<std::uint8_t>(data[1]), 5);
}

TEST(TagClipping, PileupClip_Sa_NoTrim)
{
    // sa = [3, 10, 4, 20] total = 7
    // Clip to [0, 7) = full range, no change
    const PileupClipStrategy strategy;
    TagValue val{MakeUInt8Array({3, 10, 4, 20})};

    ASSERT_TRUE(strategy.Clip(val, 0, 7, 7));
    const auto& result{std::get<TagArray>(val)};
    EXPECT_EQ(result.Count(), 4U);
    const auto data{result.Data()};
    EXPECT_EQ(static_cast<std::uint8_t>(data[0]), 3);
    EXPECT_EQ(static_cast<std::uint8_t>(data[1]), 10);
    EXPECT_EQ(static_cast<std::uint8_t>(data[2]), 4);
    EXPECT_EQ(static_cast<std::uint8_t>(data[3]), 20);
}

TEST(TagClipping, PileupClip_Sa_EmptyClip)
{
    // sa = [5, 10] total = 5
    // Clip to [2, 2) = zero length
    const PileupClipStrategy strategy;
    TagValue val{MakeUInt8Array({5, 10})};

    ASSERT_TRUE(strategy.Clip(val, 2, 0, 5));
    const auto& result{std::get<TagArray>(val)};
    EXPECT_EQ(result.Count(), 0U);
}

TEST(TagClipping, PileupClip_Sa_ManyRunsMidClip)
{
    // sa = [2, 1, 2, 2, 2, 3, 2, 4, 2, 5] total = 10
    //   pos 0-1: cov 1
    //   pos 2-3: cov 2
    //   pos 4-5: cov 3
    //   pos 6-7: cov 4
    //   pos 8-9: cov 5
    // Clip to [3, 7) = positions 3-6, clipLength=4
    //   pos 3 from run 2: length 1, cov 2
    //   pos 4-5 from run 3: length 2, cov 3
    //   pos 6 from run 4: length 1, cov 4
    // Result: [1, 2, 2, 3, 1, 4]
    const PileupClipStrategy strategy;
    TagValue val{MakeUInt8Array({2, 1, 2, 2, 2, 3, 2, 4, 2, 5})};

    ASSERT_TRUE(strategy.Clip(val, 3, 4, 10));
    const auto& result{std::get<TagArray>(val)};
    EXPECT_EQ(result.Count(), 6U);
    const auto data{result.Data()};
    EXPECT_EQ(static_cast<std::uint8_t>(data[0]), 1);
    EXPECT_EQ(static_cast<std::uint8_t>(data[1]), 2);
    EXPECT_EQ(static_cast<std::uint8_t>(data[2]), 2);
    EXPECT_EQ(static_cast<std::uint8_t>(data[3]), 3);
    EXPECT_EQ(static_cast<std::uint8_t>(data[4]), 1);
    EXPECT_EQ(static_cast<std::uint8_t>(data[5]), 4);
}

TEST(TagClipping, PileupClip_Sa_UnsupportedType_ReturnsFalse)
{
    const PileupClipStrategy strategy;
    TagValue val{std::int64_t{42}};
    EXPECT_FALSE(strategy.Clip(val, 0, 1, 1));
}

TEST(TagClipping, PileupClip_SmSx_WithPacBioDefault)
{
    // sm and sx are per-position arrays, clipped by SubstringClipStrategy via PacBioDefault
    const auto clipper{TagClipper::PacBioDefault()};

    TagMap tags;
    TagArray sm{'C'};
    for (std::uint8_t i{0}; i < 10; ++i) {
        sm.AppendUInt8(i);
    }
    tags.Set(TagKey{'s', 'm'}, std::move(sm));

    TagArray sx{'C'};
    for (std::uint8_t i{0}; i < 10; ++i) {
        sx.AppendUInt8(static_cast<std::uint8_t>(i * 2));
    }
    tags.Set(TagKey{'s', 'x'}, std::move(sx));

    clipper.ClipTags(tags, 2, 5, 10);

    const auto& smResult{std::get<TagArray>(*tags.Get(TagKey{'s', 'm'}))};
    EXPECT_EQ(smResult.Count(), 5U);
    EXPECT_EQ(static_cast<std::uint8_t>(smResult.Data()[0]), 2);

    const auto& sxResult{std::get<TagArray>(*tags.Get(TagKey{'s', 'x'}))};
    EXPECT_EQ(sxResult.Count(), 5U);
    EXPECT_EQ(static_cast<std::uint8_t>(sxResult.Data()[0]), 4);
}

TEST(TagClipping, PileupClip_SaSmSx_WithPacBioDefault)
{
    // Integration test: all three pileup tags clipped together via PacBioDefault
    const auto clipper{TagClipper::PacBioDefault()};

    TagMap tags;

    // sa: RLE coverage [5, 10, 5, 20] total = 10
    tags.Set(TagKey{'s', 'a'}, MakeUInt8Array({5, 10, 5, 20}));

    // sm: per-position matches [0,1,2,3,4,5,6,7,8,9]
    TagArray sm{'C'};
    for (std::uint8_t i{0}; i < 10; ++i) {
        sm.AppendUInt8(i);
    }
    tags.Set(TagKey{'s', 'm'}, std::move(sm));

    // sx: per-position mismatches [0,2,4,6,8,10,12,14,16,18]
    TagArray sx{'C'};
    for (std::uint8_t i{0}; i < 10; ++i) {
        sx.AppendUInt8(static_cast<std::uint8_t>(i * 2));
    }
    tags.Set(TagKey{'s', 'x'}, std::move(sx));

    // Clip to [3, 7) = positions 3-6, clipLength=4
    clipper.ClipTags(tags, 3, 4, 10);

    // sa: pos 3-4 from run 1 (len 2, cov 10), pos 5-6 from run 2 (len 2, cov 20)
    const auto& saResult{std::get<TagArray>(*tags.Get(TagKey{'s', 'a'}))};
    EXPECT_EQ(saResult.Count(), 4U);
    EXPECT_EQ(static_cast<std::uint8_t>(saResult.Data()[0]), 2);
    EXPECT_EQ(static_cast<std::uint8_t>(saResult.Data()[1]), 10);
    EXPECT_EQ(static_cast<std::uint8_t>(saResult.Data()[2]), 2);
    EXPECT_EQ(static_cast<std::uint8_t>(saResult.Data()[3]), 20);

    // sm: positions 3-6 => [3, 4, 5, 6]
    const auto& smResult{std::get<TagArray>(*tags.Get(TagKey{'s', 'm'}))};
    EXPECT_EQ(smResult.Count(), 4U);
    EXPECT_EQ(static_cast<std::uint8_t>(smResult.Data()[0]), 3);

    // sx: positions 3-6 => [6, 8, 10, 12]
    const auto& sxResult{std::get<TagArray>(*tags.Get(TagKey{'s', 'x'}))};
    EXPECT_EQ(sxResult.Count(), 4U);
    EXPECT_EQ(static_cast<std::uint8_t>(sxResult.Data()[0]), 6);
}

TEST(TagClipping, ReverseStrandMMRetainsOriginalOrientationCalls)
{
    // Reverse read: stored SEQ "TTTGGG" is revcomp of original "CCCAAA", whose three C+m
    // calls all live in the retained SEQ[3:6] window. The isReverse flag must be honoured
    // end-to-end through ClipTags or every call is wrongly dropped.
    TagMap tags;
    tags.Set(TagKey{'M', 'M'}, std::string{"C+m,0,0,0;"});

    const TagClipper clipper{TagClipper::PacBioDefault()};
    clipper.ClipTags(tags, /*clipOffset=*/3, /*clipLength=*/3, /*seqLength=*/6, "TTTGGG",
                     /*isReverse=*/true);

    const TagValue* mm{tags.Get(TagKey{'M', 'M'})};
    ASSERT_TRUE(mm);
    EXPECT_EQ(std::get<std::string>(*mm), "C+m,0,0,0;");
}

TEST(TagClipping, ClipUpdatesMnTagToClippedLength)
{
    // MN:i records the SEQ length MM/ML were produced against; htslib rejects a record
    // whose MN != l_qseq, so clipping must rewrite MN to the retained length.
    TagMap tags;
    tags.Set(TagKey{'M', 'M'}, std::string{"C+m,0;"});
    tags.Set(TagKey{'M', 'N'}, std::int64_t{10});

    const TagClipper clipper{TagClipper::PacBioDefault()};
    clipper.ClipTags(tags, /*clipOffset=*/2, /*clipLength=*/6, /*seqLength=*/10, "ACGTACGTAC");

    const TagValue* mn{tags.Get(TagKey{'M', 'N'})};
    ASSERT_TRUE(mn);
    EXPECT_EQ(std::get<std::int64_t>(*mn), 6);
}

}  // namespace Samoa
}  // namespace PacBio
