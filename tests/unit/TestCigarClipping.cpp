#include <pbsamoa/core/CigarClipping.hpp>
#include <pbsamoa/core/CigarOp.hpp>

#include <gtest/gtest.h>

namespace PacBio {
namespace Samoa {

TEST(CigarClipping, ClipToQuery_AllMatch_TrimBothEnds)
{
    // 10M, query [100,110), clip to [102,108)
    const auto cigar{ParseCigar("10M")};
    const auto result{ClipCigarToQuery(*cigar, 102, 108, 100, 110, 50, false)};
    EXPECT_EQ(CigarToString(result.cigar), "6M");
    EXPECT_EQ(result.clipOffset, 2U);
    EXPECT_EQ(result.clipLength, 6U);
    EXPECT_EQ(result.newPos, 52);
}

TEST(CigarClipping, ClipToQuery_AllMatch_NoClip)
{
    const auto cigar{ParseCigar("10M")};
    const auto result{ClipCigarToQuery(*cigar, 100, 110, 100, 110, 50, false)};
    EXPECT_EQ(CigarToString(result.cigar), "10M");
    EXPECT_EQ(result.clipOffset, 0U);
    EXPECT_EQ(result.clipLength, 10U);
    EXPECT_EQ(result.newPos, 50);
}

TEST(CigarClipping, ClipToQuery_WithInsertion)
{
    // 3M2I5M, query length = 10
    // query idx: M M M I I M M M M M
    //            0 1 2 3 4 5 6 7 8 9
    // Clip [104,110) => remove 4 from front
    // Front: eat 3M (ref+=3), eat 1 of 2I => remaining 1I + 5M
    const auto cigar{ParseCigar("3M2I5M")};
    const auto result{ClipCigarToQuery(*cigar, 104, 110, 100, 110, 50, false)};
    EXPECT_EQ(CigarToString(result.cigar), "1I5M");
    EXPECT_EQ(result.clipOffset, 4U);
    EXPECT_EQ(result.clipLength, 6U);
    EXPECT_EQ(result.newPos, 53);
}

TEST(CigarClipping, ClipToQuery_WithDeletion)
{
    // 3M2D5M, query length = 8
    // Clip [102,108) => remove 2 from front
    // Front: eat 2 of 3M (ref+=2) => 1M2D5M
    const auto cigar{ParseCigar("3M2D5M")};
    const auto result{ClipCigarToQuery(*cigar, 102, 108, 100, 108, 50, false)};
    EXPECT_EQ(CigarToString(result.cigar), "1M2D5M");
    EXPECT_EQ(result.clipOffset, 2U);
    EXPECT_EQ(result.clipLength, 6U);
    EXPECT_EQ(result.newPos, 52);
}

TEST(CigarClipping, ClipToQuery_WithSoftClips)
{
    // 2S3M1I4M2S, query=12, origQuery=[98,110)
    // Clip [100,108) => remove 2 front (2S), remove 2 back (2S)
    const auto cigar{ParseCigar("2S3M1I4M2S")};
    const auto result{ClipCigarToQuery(*cigar, 100, 108, 98, 110, 50, false)};
    EXPECT_EQ(CigarToString(result.cigar), "3M1I4M");
    EXPECT_EQ(result.clipOffset, 2U);
    EXPECT_EQ(result.clipLength, 8U);
    EXPECT_EQ(result.newPos, 50);  // soft clips don't consume reference
}

TEST(CigarClipping, ClipToQuery_Reverse)
{
    // 10M, reverse strand, clip [103,110)
    // frontRemove=3, backRemove=0
    // After reverse swap: frontRemove=0, backRemove=3
    // Result: trim 3 from back of CIGAR => 7M
    // clipOffset = 0 (nothing removed from front)
    const auto cigar{ParseCigar("10M")};
    const auto result{ClipCigarToQuery(*cigar, 103, 110, 100, 110, 50, true)};
    EXPECT_EQ(CigarToString(result.cigar), "7M");
    EXPECT_EQ(result.clipOffset, 0U);
    EXPECT_EQ(result.clipLength, 7U);
    EXPECT_EQ(result.newPos, 50);  // nothing removed from front of CIGAR
}

TEST(CigarClipping, ClipToQuery_Reverse_Symmetric)
{
    // 10M, reverse, clip [102,108) => frontRemove=2, backRemove=2
    // After swap: still frontRemove=2, backRemove=2
    const auto cigar{ParseCigar("10M")};
    const auto result{ClipCigarToQuery(*cigar, 102, 108, 100, 110, 50, true)};
    EXPECT_EQ(CigarToString(result.cigar), "6M");
    EXPECT_EQ(result.clipOffset, 2U);
    EXPECT_EQ(result.clipLength, 6U);
    EXPECT_EQ(result.newPos, 52);
}

TEST(CigarClipping, ClipToQuery_DeletionExposedAtFront)
{
    // 2D3M, query length = 3, origQuery [100,103)
    // Clip [100,102) => frontRemove=0, backRemove=1
    // Back: eat 1 of 3M => 2D2M
    const auto cigar{ParseCigar("2D3M")};
    const auto result{ClipCigarToQuery(*cigar, 100, 102, 100, 103, 50, false)};
    EXPECT_EQ(CigarToString(result.cigar), "2D2M");
    EXPECT_EQ(result.clipOffset, 0U);
    EXPECT_EQ(result.clipLength, 2U);
    EXPECT_EQ(result.newPos, 50);
}

TEST(CigarClipping, ClipToQuery_DeletionStranded)
{
    // 3M2D5M, query=8, clip front 4 => query [104,108)
    // Front: eat 3M (ref+=3), then 2D is non-query-consuming so skip (ref+=2), eat 1 of 5M (ref+=1)
    // Result: 4M
    const auto cigar{ParseCigar("3M2D5M")};
    const auto result{ClipCigarToQuery(*cigar, 104, 108, 100, 108, 50, false)};
    EXPECT_EQ(CigarToString(result.cigar), "4M");
    EXPECT_EQ(result.clipOffset, 4U);
    EXPECT_EQ(result.clipLength, 4U);
    EXPECT_EQ(result.newPos, 56);  // 3(M) + 2(D) + 1(M) = 6
}

TEST(CigarClipping, ClipToQuery_HardClips)
{
    // 3H5M3H, query=5 (hard clips don't consume query)
    // origQuery [100,105), clip [101,104)
    // Front: 3H non-query, non-ref => skip. eat 1 of 5M (ref+=1)
    // Back: 3H non-query => pop. eat 1 of 4M => 3M
    const auto cigar{ParseCigar("3H5M3H")};
    const auto result{ClipCigarToQuery(*cigar, 101, 104, 100, 105, 50, false)};
    EXPECT_EQ(CigarToString(result.cigar), "3M");
    EXPECT_EQ(result.clipOffset, 1U);
    EXPECT_EQ(result.clipLength, 3U);
    EXPECT_EQ(result.newPos, 51);
}

TEST(CigarClipping, ClipToReference_AllMatch_TrimBothEnds)
{
    // 10M at pos 50, clip to [52, 58)
    const auto cigar{ParseCigar("10M")};
    const auto result{ClipCigarToReference(*cigar, 52, 58, 50, false)};
    EXPECT_EQ(CigarToString(result.cigar), "6M");
    EXPECT_EQ(result.clipOffset, 2U);
    EXPECT_EQ(result.clipLength, 6U);
    EXPECT_EQ(result.newPos, 52);
}

TEST(CigarClipping, ClipToReference_AllMatch_NoClip)
{
    const auto cigar{ParseCigar("10M")};
    const auto result{ClipCigarToReference(*cigar, 50, 60, 50, false)};
    EXPECT_EQ(CigarToString(result.cigar), "10M");
    EXPECT_EQ(result.clipOffset, 0U);
    EXPECT_EQ(result.clipLength, 10U);
    EXPECT_EQ(result.newPos, 50);
}

TEST(CigarClipping, ClipToReference_WithInsertion_Retained)
{
    // 3M2I5M at pos 50, ref span [50,58)
    // Clip to [53, 58): remove 3 ref from front
    // Front: eat 3M (query+=3, ref-=3), insertion doesn't consume ref so stays
    // Result: 2I5M, clipOffset=3, clipLength=7
    const auto cigar{ParseCigar("3M2I5M")};
    const auto result{ClipCigarToReference(*cigar, 53, 58, 50, false)};
    EXPECT_EQ(CigarToString(result.cigar), "2I5M");
    EXPECT_EQ(result.clipOffset, 3U);
    EXPECT_EQ(result.clipLength, 7U);
    EXPECT_EQ(result.newPos, 53);
}

TEST(CigarClipping, ClipToReference_WithDeletion)
{
    // 3M2D5M at pos 50, ref span [50,60)
    // Clip to [54, 58): remove 4 ref from front, 2 from back
    // Front: eat 3M (query+=3, ref-=3), eat 1 of 2D (ref-=1) => 1D remaining
    // Result so far: 1D5M. Back: eat 2 of 5M => 1D3M
    const auto cigar{ParseCigar("3M2D5M")};
    const auto result{ClipCigarToReference(*cigar, 54, 58, 50, false)};
    EXPECT_EQ(CigarToString(result.cigar), "1D3M");
    EXPECT_EQ(result.clipOffset, 3U);
    EXPECT_EQ(result.clipLength, 3U);
    EXPECT_EQ(result.newPos, 54);
}

TEST(CigarClipping, ClipToReference_SoftClipsDropped)
{
    // 2S8M2S at pos 50 (soft clips don't affect ref span [50,58))
    // Clip to [51, 57): remove 1 ref from front, 1 from back
    // Front: 2S non-ref => pop (query+=2), eat 1 of 8M (query+=1)
    // Back: 2S non-ref => pop (query+=2), eat 1 of 7M
    // Result: 6M, clipOffset=3 (2 soft + 1 match), clipLength=6
    const auto cigar{ParseCigar("2S8M2S")};
    const auto result{ClipCigarToReference(*cigar, 51, 57, 50, false)};
    EXPECT_EQ(CigarToString(result.cigar), "6M");
    EXPECT_EQ(result.clipOffset, 3U);
    EXPECT_EQ(result.clipLength, 6U);
    EXPECT_EQ(result.newPos, 51);
}

TEST(CigarClipping, ClipToReference_ExciseFlankingInserts)
{
    // 2I3M2I at pos 50, ref span [50, 53)
    // Clip to [50, 53) with excise: no ref trimming needed, but excise flanking I's
    // Leading 2I removed (query+=2), trailing 2I removed (query+=2)
    // Result: 3M, clipOffset=2, clipLength=3
    const auto cigar{ParseCigar("2I3M2I")};
    const auto result{ClipCigarToReference(*cigar, 50, 53, 50, false, true)};
    EXPECT_EQ(CigarToString(result.cigar), "3M");
    EXPECT_EQ(result.clipOffset, 2U);
    EXPECT_EQ(result.clipLength, 3U);
    EXPECT_EQ(result.newPos, 50);
}

TEST(CigarClipping, ClipToReference_ExciseFlankingInserts_AfterRefClip)
{
    // 5M2I3M at pos 50, ref span [50,58)
    // Clip to [55, 58) with excise:
    // Front: eat 5M (query+=5, ref-=5). Now leading op is 2I => excise (query+=2)
    // Result: 3M, clipOffset=7, clipLength=3
    const auto cigar{ParseCigar("5M2I3M")};
    const auto result{ClipCigarToReference(*cigar, 55, 58, 50, false, true)};
    EXPECT_EQ(CigarToString(result.cigar), "3M");
    EXPECT_EQ(result.clipOffset, 7U);
    EXPECT_EQ(result.clipLength, 3U);
    EXPECT_EQ(result.newPos, 55);
}

TEST(CigarClipping, ClipToReference_Reverse)
{
    // 10M at pos 50, reverse strand, clip to [52, 58)
    // Front: eat 2 ref => query removed from front = 2
    // Back: eat 2 ref => query removed from back = 2
    // Reverse swap: clipOffset = queryRemovedBack = 2
    const auto cigar{ParseCigar("10M")};
    const auto result{ClipCigarToReference(*cigar, 52, 58, 50, true)};
    EXPECT_EQ(CigarToString(result.cigar), "6M");
    EXPECT_EQ(result.clipOffset, 2U);
    EXPECT_EQ(result.clipLength, 6U);
    EXPECT_EQ(result.newPos, 52);
}

TEST(CigarClipping, ClipToReference_Reverse_Asymmetric)
{
    // 10M at pos 50, reverse, clip to [53, 60) => remove 3 from front only
    // queryRemovedFront = 3, queryRemovedBack = 0
    // After reverse swap: clipOffset = 0 (back_original=0), clipLength = 7
    const auto cigar{ParseCigar("10M")};
    const auto result{ClipCigarToReference(*cigar, 53, 60, 50, true)};
    EXPECT_EQ(CigarToString(result.cigar), "7M");
    EXPECT_EQ(result.clipOffset, 0U);
    EXPECT_EQ(result.clipLength, 7U);
    EXPECT_EQ(result.newPos, 53);
}

}  // namespace Samoa
}  // namespace PacBio
