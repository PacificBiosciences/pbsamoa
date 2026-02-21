#include <pbsamoa/core/CigarOp.hpp>

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
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
    EXPECT_FALSE(CharToCigarOp('Z').has_value());
    EXPECT_FALSE(CharToCigarOp('0').has_value());
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

TEST(CigarParse, ZeroLengthOpReturnsError)
{
    EXPECT_FALSE(ParseCigar("0M").has_value());
    EXPECT_FALSE(ParseCigar("0S2M").has_value());
    EXPECT_FALSE(ParseCigar("3M0D2I").has_value());
}

TEST(CigarParse, EmptyStringReturnsEmpty)
{
    const auto cigar{ParseCigar("")};
    ASSERT_TRUE(cigar.has_value());
    EXPECT_TRUE(std::empty(*cigar));
}

TEST(CigarParse, InvalidOpCharReturnsError) { EXPECT_FALSE(ParseCigar("8Z").has_value()); }

TEST(CigarParse, MissingOpCharReturnsError) { EXPECT_FALSE(ParseCigar("8").has_value()); }

TEST(CigarParse, BadIntegerReturnsError) { EXPECT_FALSE(ParseCigar("M8").has_value()); }

}  // namespace Samoa
}  // namespace PacBio
