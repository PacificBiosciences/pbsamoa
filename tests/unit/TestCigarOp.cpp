#include <pbsamoa/core/CigarOp.hpp>

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <limits>
#include <vector>

#include <cstdint>

namespace PacBio {
namespace Samoa {

TEST(CigarOp, ConstructFromTypeAndLength)
{
    constexpr CigarOp op{CigarOpType::M, 8};
    EXPECT_EQ(op.Type(), CigarOpType::M);
    EXPECT_EQ(op.Length(), 8U);
}

TEST(CigarOp, ConstructFromRawValue)
{
    // 8M = length 8, op 0 → raw = (8 << 4) | 0 = 128
    constexpr CigarOp op{std::uint32_t{128}};
    EXPECT_EQ(op.Type(), CigarOpType::M);
    EXPECT_EQ(op.Length(), 8U);
    EXPECT_EQ(op.RawValue(), 128U);
}

TEST(CigarOp, AllOperationTypes)
{
    struct Case
    {
        CigarOpType type;
        char ch;
        std::uint8_t code;
    };

    constexpr std::array cases{
        Case{CigarOpType::M, 'M', 0}, Case{CigarOpType::I, 'I', 1},  Case{CigarOpType::D, 'D', 2},
        Case{CigarOpType::N, 'N', 3}, Case{CigarOpType::S, 'S', 4},  Case{CigarOpType::H, 'H', 5},
        Case{CigarOpType::P, 'P', 6}, Case{CigarOpType::EQ, '=', 7}, Case{CigarOpType::X, 'X', 8},
    };

    for (const auto& [type, ch, code] : cases) {
        const CigarOp op{type, 1};
        EXPECT_EQ(op.Type(), type);
        EXPECT_EQ(CigarOpToChar(type), ch) << "for code " << static_cast<int>(code);
        EXPECT_EQ(CharToCigarOp(ch), type) << "for char " << ch;
    }
}

TEST(CigarOp, CharToCigarOpRejectsInvalid)
{
    EXPECT_FALSE(CharToCigarOp('Z'));
    EXPECT_FALSE(CharToCigarOp('0'));
}

TEST(CigarOp, Equality)
{
    constexpr CigarOp a{CigarOpType::M, 8};
    constexpr CigarOp b{CigarOpType::M, 8};
    constexpr CigarOp c{CigarOpType::I, 8};
    constexpr CigarOp d{CigarOpType::M, 4};

    EXPECT_EQ(a, b);
    EXPECT_NE(a, c);
    EXPECT_NE(a, d);
}

TEST(CigarOp, ConsumesQueryAndReference)
{
    // M,=,X consume both; I,S consume query only; D,N consume ref only; H,P consume neither
    EXPECT_TRUE(ConsumesQuery(CigarOpType::M));
    EXPECT_TRUE(ConsumesReference(CigarOpType::M));
    EXPECT_TRUE(ConsumesQuery(CigarOpType::I));
    EXPECT_FALSE(ConsumesReference(CigarOpType::I));
    EXPECT_FALSE(ConsumesQuery(CigarOpType::D));
    EXPECT_TRUE(ConsumesReference(CigarOpType::D));
    EXPECT_FALSE(ConsumesQuery(CigarOpType::H));
    EXPECT_FALSE(ConsumesReference(CigarOpType::H));
}

TEST(CigarParse, SpecExample)
{
    // r001: 8M2I4M1D3M
    const std::vector<CigarOp> cigar{*ParseCigar("8M2I4M1D3M")};
    ASSERT_EQ(std::size(cigar), 5U);
    EXPECT_EQ(cigar[0], CigarOp(CigarOpType::M, 8));
    EXPECT_EQ(cigar[1], CigarOp(CigarOpType::I, 2));
    EXPECT_EQ(cigar[2], CigarOp(CigarOpType::M, 4));
    EXPECT_EQ(cigar[3], CigarOp(CigarOpType::D, 1));
    EXPECT_EQ(cigar[4], CigarOp(CigarOpType::M, 3));
}

TEST(CigarParse, StarReturnsEmpty)
{
    const std::vector<CigarOp> cigar{*ParseCigar("*")};
    EXPECT_TRUE(std::empty(cigar));
}

TEST(CigarParse, SingleOp)
{
    const std::vector<CigarOp> cigar{*ParseCigar("100M")};
    ASSERT_EQ(std::size(cigar), 1U);
    EXPECT_EQ(cigar[0], CigarOp(CigarOpType::M, 100));
}

TEST(CigarParse, AllOpTypes)
{
    const std::vector<CigarOp> cigar{*ParseCigar("1M1I1D1N1S1H1P1=1X")};
    ASSERT_EQ(std::size(cigar), 9U);
    EXPECT_EQ(cigar[0].Type(), CigarOpType::M);
    EXPECT_EQ(cigar[8].Type(), CigarOpType::X);
}

TEST(CigarParse, RejectsOpLengthExceeding28BitField)
{
    // The BAM CIGAR field is len<<4|op; a length >= 2^28 would clobber the op code.
    EXPECT_TRUE(ParseCigar("268435455M").has_value());   // 2^28 - 1: largest valid
    EXPECT_FALSE(ParseCigar("268435456M").has_value());  // 2^28: out of range
}

TEST(CigarToString, RoundTrip)
{
    const std::string original{"8M2I4M1D3M"};
    const std::vector<CigarOp> cigar{*ParseCigar(original)};
    EXPECT_EQ(CigarToString(cigar), original);
}

TEST(CigarToString, EmptyReturnsStar) { EXPECT_EQ(CigarToString({}), "*"); }

TEST(CigarLength, SpecExampleReferenceLength)
{
    // 8M2I4M1D3M: ref = 8+4+1+3 = 16
    const std::vector<CigarOp> cigar{*ParseCigar("8M2I4M1D3M")};
    EXPECT_EQ(ReferenceLength(cigar), 16);
}

TEST(CigarLength, SpecExampleQueryLength)
{
    // 8M2I4M1D3M: query = 8+2+4+3 = 17
    const std::vector<CigarOp> cigar{*ParseCigar("8M2I4M1D3M")};
    EXPECT_EQ(QueryLength(cigar), 17);
}

TEST(CigarLength, WithSoftClip)
{
    // 3S6M1P1I4M: query = 3+6+1+4 = 14, ref = 6+4 = 10
    const std::vector<CigarOp> cigar{*ParseCigar("3S6M1P1I4M")};
    EXPECT_EQ(QueryLength(cigar), 14);
    EXPECT_EQ(ReferenceLength(cigar), 10);
}

TEST(Reg2Bin, SmallRegion)
{
    // [0, 100) fits in a single 16kbp bin → level 5
    EXPECT_EQ(Reg2Bin(0, 100), 4681U);
}

TEST(Reg2Bin, CrossesBinBoundary)
{
    // [0, 16384) spans exactly one 16kbp bin
    EXPECT_EQ(Reg2Bin(0, 16384), 4681U);
    // [0, 16385) crosses into second bin → goes up a level
    EXPECT_EQ(Reg2Bin(0, 16385), 585U);
}

TEST(Reg2Bin, LargeRegion)
{
    // [0, 2^29) — spans everything → bin 0
    EXPECT_EQ(Reg2Bin(0, 1 << 29), 0U);
}

TEST(Reg2Bins, AlwaysContainsBinZero)
{
    const auto bins = Reg2Bins(0, 100);
    EXPECT_FALSE(std::empty(bins));
    EXPECT_EQ(bins[0], 0U);
}

TEST(Reg2Bins, IncludesReg2BinResult)
{
    const std::int32_t beg{100};
    const std::int32_t end{200};
    const std::uint16_t expected{Reg2Bin(beg, end)};
    const auto bins = Reg2Bins(beg, end);
    EXPECT_NE(std::ranges::find(bins, expected), std::ranges::end(bins));
}

TEST(Reg2Bins, LargeRegionProducesManyBins)
{
    // A region spanning several 16kbp windows should yield many bins
    const auto bins = Reg2Bins(0, 1 << 20);
    EXPECT_GT(std::size(bins), 10U);
}

TEST(CigarParse, ZeroLengthOpAccepted)
{
    // Zero-length ops are grammar-valid and accepted by htslib; pbsamoa matches that.
    EXPECT_TRUE(ParseCigar("0M").has_value());
    EXPECT_TRUE(ParseCigar("0S2M").has_value());
    const auto cigar{ParseCigar("3M0D2I")};
    ASSERT_TRUE(cigar);
    ASSERT_EQ(std::size(*cigar), 3U);
    EXPECT_EQ((*cigar)[1], CigarOp(CigarOpType::D, 0));
}

TEST(CigarParse, EmptyStringReturnsEmpty)
{
    const auto cigar{ParseCigar("")};
    ASSERT_TRUE(cigar);
    EXPECT_TRUE(std::empty(*cigar));
}

TEST(CigarParse, InvalidOpCharReturnsError) { EXPECT_FALSE(ParseCigar("8Z")); }

TEST(CigarParse, MissingOpCharReturnsError) { EXPECT_FALSE(ParseCigar("8")); }

TEST(CigarParse, BadIntegerReturnsError) { EXPECT_FALSE(ParseCigar("M8")); }

// Reg2Bins must clamp 'end' to 1<<29 so queries beyond the max binnable
// coordinate return the same bin set as end==1<<29, matching htslib reg2bins.
TEST(Reg2Bins, EndBeyondMaxCoordClampsTo1Shl29)
{
    constexpr std::int32_t MAX_COORD{1 << 29};
    const auto binsAtMax{Reg2Bins(0, MAX_COORD)};
    const auto binsOverflow{Reg2Bins(0, MAX_COORD + 1)};
    const auto binsWayOver{Reg2Bins(0, std::numeric_limits<std::int32_t>::max())};
    EXPECT_EQ(binsAtMax, binsOverflow);
    EXPECT_EQ(binsAtMax, binsWayOver);
}

// Reg2Bins must return empty when beg >= end (after clamping), matching
// htslib's early-return guard.
TEST(Reg2Bins, BegAtOrAfterClampedEndReturnsEmpty)
{
    constexpr std::int32_t MAX_COORD{1 << 29};
    // beg == end: half-open interval is empty
    EXPECT_TRUE(std::empty(Reg2Bins(100, 100)));
    // beg > end after clamp: e.g. beg == 1<<29, end == 1<<29 (clamped)
    EXPECT_TRUE(std::empty(Reg2Bins(MAX_COORD, MAX_COORD + 1)));
}

// ComputeCigarStats must match pbcopper CigarOpsCalculator definitions exactly,
// so an aligner emitting de/mc/identity tags via pbsamoa matches pbmm2 output.
TEST(CigarStats, MatchesPbcopperCalculatorDefinitions)
{
    // 50= 2X 3I 50= 4D 50=  ->  match=150, mismatch=2, ins=3 (1 event), del=4 (1 event)
    const auto cigar{*ParseCigar("50=2X3I50=4D50=")};
    const CigarStats stats{ComputeCigarStats(cigar)};
    EXPECT_EQ(stats.Mismatches, 2);
    EXPECT_EQ(stats.Insertions, 3);
    EXPECT_EQ(stats.Deletions, 4);
    EXPECT_EQ(stats.InsertionEvents, 1);
    EXPECT_EQ(stats.DeletionEvents, 1);
    // NumAlignedBases = match + insertions + mismatch = 150 + 3 + 2
    EXPECT_EQ(stats.NumAlignedBases, 155);
    EXPECT_EQ(stats.Span, 155);
    // Identity = 100 * match / (match + mismatch + del + ins)
    EXPECT_NEAR(stats.Identity, 100.0 * 150.0 / 159.0, 1e-9);
    // IdentityGapComp = 100 * match / (match + mismatch + delEvents + insEvents)
    EXPECT_NEAR(stats.IdentityGapComp, 100.0 * 150.0 / 154.0, 1e-9);
}

TEST(CigarStats, TreatsMOpAsMatch)
{
    const auto cigar{*ParseCigar("10M")};
    const CigarStats stats{ComputeCigarStats(cigar)};
    EXPECT_EQ(stats.NumAlignedBases, 10);
    EXPECT_DOUBLE_EQ(stats.Identity, 100.0);
}

TEST(CigarStats, ExcludesClipsAndRefSkips)
{
    const auto cigar{*ParseCigar("5S10=3N10=5H")};
    const CigarStats stats{ComputeCigarStats(cigar)};
    EXPECT_EQ(stats.NumAlignedBases, 20);
    EXPECT_DOUBLE_EQ(stats.Identity, 100.0);
}

TEST(CigarStats, EmptyCigarYieldsZeroIdentity)
{
    const CigarStats stats{ComputeCigarStats(CigarView{})};
    EXPECT_EQ(stats.NumAlignedBases, 0);
    EXPECT_DOUBLE_EQ(stats.Identity, 0.0);
    EXPECT_DOUBLE_EQ(stats.IdentityGapComp, 0.0);
}

TEST(CountMismatches, SumsXOps)
{
    const auto cigar{*ParseCigar("10=2X5=3X")};
    EXPECT_EQ(CountMismatches(cigar), 5);
}

}  // namespace Samoa
}  // namespace PacBio
