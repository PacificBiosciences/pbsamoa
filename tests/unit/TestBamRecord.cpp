#include <pbsamoa/core/BamRecord.hpp>
#include <pbsamoa/core/Endian.hpp>
#include <pbsamoa/core/RawRecord.hpp>
#include <pbsamoa/core/Sequence.hpp>
#include <pbsamoa/core/Tags.hpp>

#include <gtest/gtest.h>

#include <algorithm>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {

TEST(BamRecord, DefaultConstruction)
{
    const BamRecord rec;
    EXPECT_TRUE(std::empty(rec.Name()));
    EXPECT_EQ(rec.Flag(), 0U);
    EXPECT_EQ(rec.RefId(), -1);
    EXPECT_EQ(rec.Pos(), -1);
    EXPECT_EQ(rec.MapQ(), 0U);
    EXPECT_TRUE(std::empty(rec.Cigar()));
    EXPECT_EQ(rec.NextRefId(), -1);
    EXPECT_EQ(rec.NextPos(), -1);
    EXPECT_EQ(rec.Tlen(), 0);
    EXPECT_TRUE(std::empty(rec.Sequence()));
    EXPECT_TRUE(std::empty(rec.Qualities()));
    EXPECT_TRUE(rec.Tags().Empty());
}

TEST(BamRecord, SettersAndGetters)
{
    BamRecord rec;
    rec.Name("r001")
        .Flag(99)
        .RefId(0)
        .Pos(6)
        .MapQ(30)
        .Cigar(*ParseCigar("8M2I4M1D3M"))
        .NextRefId(0)
        .NextPos(36)
        .Tlen(39)
        .Sequence("TTAGATAAAGGATACTG")
        .Qualities({30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30});

    EXPECT_EQ(rec.Name(), "r001");
    EXPECT_EQ(rec.Flag(), 99U);
    EXPECT_EQ(rec.RefId(), 0);
    EXPECT_EQ(rec.Pos(), 6);
    EXPECT_EQ(rec.MapQ(), 30U);
    EXPECT_EQ(std::size(rec.Cigar()), 5U);
    EXPECT_EQ(rec.NextRefId(), 0);
    EXPECT_EQ(rec.NextPos(), 36);
    EXPECT_EQ(rec.Tlen(), 39);
    EXPECT_EQ(rec.Sequence(), "TTAGATAAAGGATACTG");
    EXPECT_EQ(std::size(rec.Qualities()), 17U);
}

TEST(BamRecord, TagAccess)
{
    BamRecord rec;
    TagMap tags;
    tags.Set(TagKey{'N', 'M'}, TagValue{std::int64_t{1}});
    rec.Tags(std::move(tags));

    EXPECT_TRUE(rec.Tags().Contains(TagKey{'N', 'M'}));
    EXPECT_EQ(std::get<std::int64_t>(*rec.Tags().Get(TagKey{'N', 'M'})), 1);
}

TEST(BamRecord, IsMapped)
{
    BamRecord mapped;
    mapped.Flag(0);
    EXPECT_TRUE(mapped.IsMapped());

    BamRecord unmapped;
    unmapped.Flag(4);  // 0x4 = segment unmapped
    EXPECT_FALSE(unmapped.IsMapped());
}

TEST(BamRecord, IsReverse)
{
    BamRecord fwd;
    fwd.Flag(0);
    EXPECT_FALSE(fwd.IsReverseStrand());

    BamRecord rev;
    rev.Flag(16);  // 0x10 = reverse strand
    EXPECT_TRUE(rev.IsReverseStrand());
}

TEST(BamRecord, ReferenceEnd)
{
    BamRecord rec;
    rec.Pos(6).Cigar(*ParseCigar("8M2I4M1D3M"));
    // RefLen = 8+4+1+3 = 16, end = 6 + 16 = 22
    EXPECT_EQ(rec.ReferenceEnd(), 22);
}

TEST(BamRecord, SerializeToBamBasic)
{
    BamRecord rec;
    rec.Name("r001")
        .Flag(99)
        .RefId(0)
        .Pos(6)
        .MapQ(30)
        .Cigar(*ParseCigar("4M"))
        .NextRefId(0)
        .NextPos(36)
        .Tlen(39)
        .Sequence("ACGT")
        .Qualities({30, 30, 30, 30});

    const std::vector<std::byte> bam{rec.SerializeToBam()};
    EXPECT_FALSE(std::empty(bam));

    // Verify fixed fields at known offsets
    EXPECT_EQ(ReadI32LE(std::data(bam) + 0), 0);
    EXPECT_EQ(ReadI32LE(std::data(bam) + 4), 6);

    const std::uint8_t nameLen{static_cast<std::uint8_t>(bam[8])};
    EXPECT_EQ(nameLen, 5U);  // "r001" + NUL

    EXPECT_EQ(static_cast<std::uint8_t>(bam[9]), 30U);  // mapq

    EXPECT_EQ(ReadU16LE(std::data(bam) + 12), 1U);
    EXPECT_EQ(ReadU16LE(std::data(bam) + 14), 99U);
    EXPECT_EQ(ReadU32LE(std::data(bam) + 16), 4U);
}

TEST(BamRecord, SerializeToBamRoundTrip)
{
    BamRecord rec;
    rec.Name("read1")
        .Flag(0)
        .RefId(0)
        .Pos(9)
        .MapQ(40)
        .Cigar(*ParseCigar("4M"))
        .Sequence("ACGT")
        .Qualities({30, 30, 30, 30});

    TagMap tags;
    tags.Set(TagKey{'N', 'M'}, TagValue{std::int64_t{1}});
    rec.Tags(std::move(tags));

    const std::vector<std::byte> bam{rec.SerializeToBam()};

    const std::uint8_t nameLen{static_cast<std::uint8_t>(bam[8])};
    const std::string name{reinterpret_cast<const char*>(std::data(bam) + 32),
                           static_cast<std::size_t>(nameLen - 1)};
    EXPECT_EQ(name, "read1");

    const std::uint32_t lSeq{4};
    const std::size_t seqOffset{32U + nameLen + 4U};  // past name + 1 cigar op
    const std::string seq{UnpackSequence(std::span{bam}.subspan(seqOffset, 2), lSeq)};
    EXPECT_EQ(seq, "ACGT");
}

TEST(BamRecord, SerializeUnmappedRecord)
{
    BamRecord rec;
    rec.Name("unmapped").Flag(4).RefId(-1).Pos(-1);

    const std::vector<std::byte> bam{rec.SerializeToBam()};
    EXPECT_FALSE(std::empty(bam));

    EXPECT_EQ(ReadI32LE(std::data(bam) + 0), -1);
}

TEST(BamRecord, SerializeNoSequence)
{
    BamRecord rec;
    rec.Name("noSeq").Flag(4).RefId(-1).Pos(-1);

    const std::vector<std::byte> bam{rec.SerializeToBam()};
    EXPECT_EQ(ReadU32LE(std::data(bam) + 16), 0U);
}

TEST(BamRecord, SerializeToBamUsesLittleEndianFixedFields)
{
    BamRecord rec;
    rec.Name("le")
        .Flag(0x1234)
        .RefId(0x01020304)
        .Pos(0x0A0B0C0D)
        .MapQ(0x1E)
        .Cigar(*ParseCigar("1M"))
        .NextRefId(0x11121314)
        .NextPos(0x21222324)
        .Tlen(0x31323334)
        .Sequence("A")
        .Qualities({30});

    const std::vector<std::byte> bam{rec.SerializeToBam()};
    ASSERT_GE(std::size(bam), 32u);

    // refId
    EXPECT_EQ(bam[0], static_cast<std::byte>(0x04));
    EXPECT_EQ(bam[1], static_cast<std::byte>(0x03));
    EXPECT_EQ(bam[2], static_cast<std::byte>(0x02));
    EXPECT_EQ(bam[3], static_cast<std::byte>(0x01));
    // pos
    EXPECT_EQ(bam[4], static_cast<std::byte>(0x0D));
    EXPECT_EQ(bam[5], static_cast<std::byte>(0x0C));
    EXPECT_EQ(bam[6], static_cast<std::byte>(0x0B));
    EXPECT_EQ(bam[7], static_cast<std::byte>(0x0A));
    // flag
    EXPECT_EQ(bam[14], static_cast<std::byte>(0x34));
    EXPECT_EQ(bam[15], static_cast<std::byte>(0x12));
}

TEST(BamRecord, NameLengthValidation)
{
    BamRecord rec;

    // Max allowed: 254 characters
    const std::string maxName(254, 'A');
    EXPECT_NO_THROW(rec.Name(maxName));
    EXPECT_EQ(std::size(rec.Name()), 254u);

    // Over limit: 255 characters
    const std::string tooLong(255, 'A');
    EXPECT_THROW(rec.Name(tooLong), std::invalid_argument);
}

TEST(BamRecord, SerializeToBamNameOverflowThrows)
{
    BamRecord rec;
    // Construct with short name, then verify serialize works
    rec.Name("ok").Flag(0).Sequence("A");
    EXPECT_NO_THROW(static_cast<void>(rec.SerializeToBam()));
}

TEST(BamRecord, EmptyNameSerializesAsStar)
{
    // An empty QNAME is written as the placeholder "*" (htslib bam_set1).
    BamRecord rec;
    rec.Name("").Flag(4).Sequence("");
    const std::vector<std::byte> bytes{rec.SerializeToBam()};

    const std::uint8_t lReadName{static_cast<std::uint8_t>(bytes[8])};
    ASSERT_EQ(lReadName, 2U);  // '*' + NUL
    EXPECT_EQ(static_cast<char>(bytes[32]), '*');
    EXPECT_EQ(static_cast<char>(bytes[33]), '\0');
}

TEST(BamRecord, LongCigarRoundTripsViaCgTag)
{
    // A CIGAR with > 65535 ops must serialize as the kSmN placeholder + CG:B,I tag and
    // decode back to the full CIGAR (SAMv1 §4.2.2).
    constexpr std::uint32_t opCount{70000};
    std::vector<CigarOp> cigar;
    cigar.reserve(opCount);
    for (std::uint32_t i{0}; i < opCount; ++i) {
        cigar.emplace_back(CigarOpType::M, 1);  // each 1M consumes one query + ref base
    }

    BamRecord rec;
    rec.Name("longcig")
        .Flag(0)
        .RefId(0)
        .Pos(0)
        .Cigar(std::move(cigar))
        .Sequence(std::string(opCount, 'A'));

    const std::vector<std::byte> bytes{rec.SerializeToBam()};

    // On-disk n_cigar_op must be the 2-op placeholder, not the real count.
    EXPECT_EQ(ReadU16LE(std::data(bytes) + 12), 2U);

    // Decode: RawRecord must expand the CG tag back to the full CIGAR and drop the CG tag.
    const RawRecord decoded{std::span<const std::byte>{bytes}};
    EXPECT_EQ(std::size(decoded.CigarOps()), opCount);
    EXPECT_FALSE(decoded.ParseTags().Contains(TagKey{'C', 'G'}));
}

}  // namespace Samoa
}  // namespace PacBio
