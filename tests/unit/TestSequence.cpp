#include <pbsamoa/core/Sequence.hpp>

#include <gtest/gtest.h>

#include <string>
#include <vector>

#include <cstddef>

namespace PacBio {
namespace Samoa {

TEST(SequenceEncoding, PackUnpackRoundTrip)
{
    const std::string seq{"ACGTACGT"};
    const std::vector<std::byte> packed{PackSequence(seq)};
    EXPECT_EQ(std::size(packed), 4U);  // 8 bases / 2 = 4 bytes
    EXPECT_EQ(UnpackSequence(packed, 8), seq);
}

TEST(SequenceEncoding, OddLength)
{
    const std::string seq{"ACG"};
    const std::vector<std::byte> packed{PackSequence(seq)};
    EXPECT_EQ(std::size(packed), 2U);  // (3+1)/2 = 2 bytes
    EXPECT_EQ(UnpackSequence(packed, 3), seq);
}

TEST(SequenceEncoding, SingleBase)
{
    const std::string seq{"T"};
    const std::vector<std::byte> packed{PackSequence(seq)};
    EXPECT_EQ(std::size(packed), 1U);
    EXPECT_EQ(UnpackSequence(packed, 1), seq);
}

TEST(SequenceEncoding, Empty)
{
    const std::vector<std::byte> packed{PackSequence("")};
    EXPECT_TRUE(std::empty(packed));
    EXPECT_EQ(UnpackSequence({}, 0), "");
}

TEST(SequenceEncoding, AllAmbiguityCodes)
{
    const std::string seq{"=ACMGRSVTWYHKDBN"};
    const std::vector<std::byte> packed{PackSequence(seq)};
    EXPECT_EQ(std::size(packed), 8U);  // 16 bases / 2
    EXPECT_EQ(UnpackSequence(packed, 16), seq);
}

TEST(SequenceEncoding, SpecExampleSequence)
{
    // r001: TTAGATAAAGGATACTG (17 bases)
    const std::string seq{"TTAGATAAAGGATACTG"};
    const std::vector<std::byte> packed{PackSequence(seq)};
    EXPECT_EQ(std::size(packed), 9U);  // (17+1)/2
    EXPECT_EQ(UnpackSequence(packed, 17), seq);
}

TEST(SequenceView, AccessByIndex)
{
    const std::string seq{"ACGT"};
    const std::vector<std::byte> packed{PackSequence(seq)};
    const SequenceView view{packed, 4};

    EXPECT_EQ(view.Size(), 4U);
    EXPECT_EQ(view[0], 'A');
    EXPECT_EQ(view[1], 'C');
    EXPECT_EQ(view[2], 'G');
    EXPECT_EQ(view[3], 'T');
}

TEST(SequenceView, ToString)
{
    const std::string seq{"TTAGATAAAGGATACTG"};
    const std::vector<std::byte> packed{PackSequence(seq)};
    const SequenceView view{packed, 17};
    EXPECT_EQ(view.ToString(), seq);
}

TEST(SequenceEncoding, PackSequenceIntoRoundTrip)
{
    const std::string seq{"ACGTACGT"};
    const std::size_t packedLen{(std::size(seq) + 1) / 2};
    std::vector<std::byte> buf(packedLen);
    PackSequenceInto(seq, std::data(buf));
    EXPECT_EQ(UnpackSequence(buf, static_cast<std::uint32_t>(std::size(seq))), seq);
}

TEST(SequenceEncoding, PackSequenceIntoOddLength)
{
    const std::string seq{"ACG"};
    std::vector<std::byte> buf(2);
    PackSequenceInto(seq, std::data(buf));
    EXPECT_EQ(UnpackSequence(buf, 3), seq);
}

TEST(SequenceEncoding, PackSequenceIntoSingleBase)
{
    const std::string seq{"T"};
    std::vector<std::byte> buf(1);
    PackSequenceInto(seq, std::data(buf));
    EXPECT_EQ(UnpackSequence(buf, 1), seq);
}

TEST(Sequence, ReverseComplement_basic)
{
    EXPECT_EQ(ReverseComplement("ACGT"), "ACGT");
    EXPECT_EQ(ReverseComplement("AAAA"), "TTTT");
    EXPECT_EQ(ReverseComplement("A"), "T");
    EXPECT_EQ(ReverseComplement(""), "");
    EXPECT_EQ(ReverseComplement("GATTACA"), "TGTAATC");
}

TEST(Sequence, ReverseComplement_preserves_N)
{
    EXPECT_EQ(ReverseComplement("ACNGT"), "ACNGT");
    EXPECT_EQ(ReverseComplement("NNN"), "NNN");
}

TEST(Sequence, ReverseComplement_preserves_equals_base)
{
    // '=' (match-to-reference) is self-complementary and must round-trip, not collapse to N.
    EXPECT_EQ(ReverseComplement("A=C"), "G=T");
}

TEST(Sequence, ReverseComplement_lowercase_iupac)
{
    // Lowercase IUPAC ambiguity codes complement case-preservingly (not collapse to N).
    EXPECT_EQ(ReverseComplement("rymkswbvdh"), "dhbvwsmkry");
}

TEST(Sequence, ReverseComplementInPlace_basic)
{
    std::string seq{"ACGT"};
    ReverseComplementInPlace(seq);
    EXPECT_EQ(seq, "ACGT");

    std::string seq2{"GATTACA"};
    ReverseComplementInPlace(seq2);
    EXPECT_EQ(seq2, "TGTAATC");
}

}  // namespace Samoa
}  // namespace PacBio
