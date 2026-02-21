#include <pbsamoa/core/RawRecord.hpp>

#include <pbsamoa/core/BamRecord.hpp>
#include <pbsamoa/core/CigarOp.hpp>
#include <pbsamoa/core/Tags.hpp>

#include <gtest/gtest.h>

#include <span>
#include <vector>

#include <cstddef>

namespace PacBio {
namespace Samoa {

namespace {

std::vector<std::byte> MakeTestRecordBytes()
{
    BamRecord rec;
    rec.Name("r001")
        .Flag(99)
        .RefId(0)
        .Pos(6)
        .MapQ(30)
        .Cigar(ParseCigar("8M2I4M1D3M"))
        .NextRefId(0)
        .NextPos(36)
        .Tlen(39)
        .Sequence("TTAGATAAAGGATACTG")
        .Qualities({30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30});

    TagMap tags;
    tags.Set(TagKey{'N', 'M'}, TagValue{std::int64_t{1}});
    rec.Tags(std::move(tags));

    return rec.SerializeToBam();
}

}  // namespace

TEST(RawRecord, ConstructFromSpan)
{
    const std::vector<std::byte> data{MakeTestRecordBytes()};
    const RawRecord view{std::span<const std::byte>{data}};
    EXPECT_FALSE(view.RawData().empty());
}

TEST(RawRecord, OwnsData)
{
    std::vector<std::byte> data{MakeTestRecordBytes()};
    const RawRecord view{std::span<const std::byte>{data}};

    // Overwrite original data — view must still work because it owns a copy
    std::ranges::fill(data, std::byte{0xFF});

    EXPECT_EQ(view.RefId(), 0);
    EXPECT_EQ(view.Pos(), 6);
    EXPECT_EQ(view.Name(), "r001");
}

TEST(RawRecord, FieldAccessors)
{
    const std::vector<std::byte> data{MakeTestRecordBytes()};
    const RawRecord view{std::span<const std::byte>{data}};

    EXPECT_EQ(view.RefId(), 0);
    EXPECT_EQ(view.Pos(), 6);
    EXPECT_EQ(view.MapQ(), 30u);
    EXPECT_EQ(view.Flag(), 99u);
    EXPECT_EQ(view.NextRefId(), 0);
    EXPECT_EQ(view.NextPos(), 36);
    EXPECT_EQ(view.Tlen(), 39);
    EXPECT_EQ(view.Name(), "r001");
    EXPECT_EQ(view.SeqLength(), 17);
    EXPECT_EQ(view.CigarOpCount(), 5u);
}

TEST(RawRecord, ToOwnedProducesEquivalentRecord)
{
    const std::vector<std::byte> data{MakeTestRecordBytes()};
    const RawRecord view{std::span<const std::byte>{data}};

    const BamRecord owned{view.ToOwned()};

    EXPECT_EQ(owned.Name(), "r001");
    EXPECT_EQ(owned.Flag(), 99u);
    EXPECT_EQ(owned.Pos(), 6);
    EXPECT_EQ(owned.Tlen(), 39);
    EXPECT_EQ(owned.Sequence(), "TTAGATAAAGGATACTG");
}

TEST(RawRecord, MoveSemantics)
{
    const std::vector<std::byte> data{MakeTestRecordBytes()};
    RawRecord view1{std::span<const std::byte>{data}};

    const RawRecord view2{std::move(view1)};
    EXPECT_EQ(view2.Pos(), 6);
    EXPECT_EQ(view2.Name(), "r001");
}

TEST(RawRecord, TooSmallForFixedFieldsThrows)
{
    const std::vector<std::byte> data(16, std::byte{0});
    EXPECT_THROW(RawRecord(std::span<const std::byte>{data}), std::invalid_argument);
}

TEST(RawRecord, EmptyDataThrows)
{
    const std::vector<std::byte> data;
    EXPECT_THROW(RawRecord(std::span<const std::byte>{data}), std::invalid_argument);
}

TEST(RawRecord, ZeroNameLengthThrows)
{
    // Create 32 bytes (minimum fixed fields) with l_read_name = 0 at offset 8
    std::vector<std::byte> data(32, std::byte{0});
    EXPECT_THROW(RawRecord(std::span<const std::byte>{data}), std::invalid_argument);
}

TEST(RawRecord, TruncatedVariableLengthFieldsThrows)
{
    // Create valid fixed fields but claim name length > remaining buffer
    std::vector<std::byte> data(33, std::byte{0});
    data[8] = std::byte{10};  // l_read_name = 10, but only 1 byte after fixed fields
    EXPECT_THROW(RawRecord(std::span<const std::byte>{data}), std::invalid_argument);
}

TEST(RawRecord, ExactMinimumSizeIsValid)
{
    // 32 fixed bytes + 1 byte name (l_read_name=1) + no CIGAR + no seq + no qual
    std::vector<std::byte> data(33, std::byte{0});
    data[8] = std::byte{1};  // l_read_name = 1 (just the NUL)
    EXPECT_NO_THROW(RawRecord(std::span<const std::byte>{data}));
}

}  // namespace Samoa
}  // namespace PacBio
