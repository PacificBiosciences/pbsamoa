#include <pbsamoa/core/RawRecord.hpp>

#include <pbsamoa/core/BamRecord.hpp>
#include <pbsamoa/core/CigarOp.hpp>

#include <gtest/gtest.h>

#include <algorithm>
#include <span>
#include <vector>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {

namespace {

std::uint32_t ReadU32(const std::byte* data)
{
    std::uint32_t value{0};
    std::ranges::copy_n(data, sizeof(value), reinterpret_cast<std::byte*>(&value));
    return value;
}

void WriteU32(std::byte* data, std::uint32_t value)
{
    std::ranges::copy_n(reinterpret_cast<const std::byte*>(&value), sizeof(value), data);
}

// Build a buffer containing multiple serialized records (with block_size
// prefixes)
std::vector<std::byte> BuildMultiRecordBuffer(std::vector<std::uint32_t>& recordOffsets)
{
    std::vector<std::byte> buffer;

    BamRecord rec1;
    rec1.Name("read1")
        .Flag(0)
        .RefId(0)
        .Pos(100)
        .MapQ(30)
        .Cigar(*ParseCigar("10M"))
        .Sequence("ACGTACGTAC")
        .Qualities({30, 30, 30, 30, 30, 30, 30, 30, 30, 30});

    BamRecord rec2;
    rec2.Name("read2")
        .Flag(16)
        .RefId(0)
        .Pos(200)
        .MapQ(40)
        .Cigar(*ParseCigar("5M"))
        .Sequence("TTTTT")
        .Qualities({40, 40, 40, 40, 40});

    BamRecord rec3;
    rec3.Name("read3").Flag(4).RefId(-1).Pos(-1).Sequence("GGG");

    for (const BamRecord& rec : {rec1, rec2, rec3}) {
        const std::vector<std::byte> recBytes{rec.SerializeToBam()};
        const std::uint32_t blockSize = std::size(recBytes);

        // Record offset is position of the record data (after block_size)
        recordOffsets.push_back(std::size(buffer) + 4);

        std::vector<std::byte> tmp(4);
        WriteU32(std::data(tmp), blockSize);
        buffer.insert(std::ranges::end(buffer), std::ranges::begin(tmp), std::ranges::end(tmp));
        buffer.insert(std::ranges::end(buffer), std::ranges::begin(recBytes),
                      std::ranges::end(recBytes));
    }

    return buffer;
}

}  // namespace

TEST(RawRecordBatch, ConstructAndIterate)
{
    std::vector<std::uint32_t> offsets;
    std::vector<std::byte> buffer{BuildMultiRecordBuffer(offsets)};

    // Build record sizes from offsets
    std::vector<RawRecordBatch::RecordExtent> extents;
    for (std::size_t i{0}; i < std::size(offsets); ++i) {
        const std::uint32_t dataOffset{offsets[i]};
        const std::uint32_t blockSize{ReadU32(std::data(buffer) + dataOffset - 4)};
        extents.push_back({dataOffset, blockSize});
    }

    const RawRecordBatch batch{std::move(buffer), std::move(extents)};

    ASSERT_EQ(batch.RecordCount(), 3U);
    const RawRecord view0{batch.RecordData(0)};
    const RawRecord view1{batch.RecordData(1)};
    const RawRecord view2{batch.RecordData(2)};
    EXPECT_EQ(view0.Name(), "read1");
    EXPECT_EQ(view1.Name(), "read2");
    EXPECT_EQ(view2.Name(), "read3");
}

TEST(RawRecordBatch, RecordFieldAccess)
{
    std::vector<std::uint32_t> offsets;
    std::vector<std::byte> buffer{BuildMultiRecordBuffer(offsets)};

    std::vector<RawRecordBatch::RecordExtent> extents;
    for (std::size_t i{0}; i < std::size(offsets); ++i) {
        const std::uint32_t dataOffset{offsets[i]};
        const std::uint32_t blockSize{ReadU32(std::data(buffer) + dataOffset - 4)};
        extents.push_back({dataOffset, blockSize});
    }

    const RawRecordBatch batch{std::move(buffer), std::move(extents)};

    const RawRecord view0{batch.RecordData(0)};
    const RawRecord view1{batch.RecordData(1)};
    const RawRecord view2{batch.RecordData(2)};

    EXPECT_EQ(view0.Pos(), 100);
    EXPECT_TRUE(view0.IsMapped());

    EXPECT_EQ(view1.Pos(), 200);
    EXPECT_TRUE(view1.IsReverseStrand());

    EXPECT_EQ(view2.RefId(), -1);
    EXPECT_FALSE(view2.IsMapped());
}

TEST(RawRecordBatch, EmptyBatch)
{
    const RawRecordBatch batch{std::vector<std::byte>{}, {}};
    EXPECT_EQ(batch.RecordCount(), 0U);
}

TEST(RawRecordBatch, RecordCount)
{
    std::vector<std::uint32_t> offsets;
    std::vector<std::byte> buffer{BuildMultiRecordBuffer(offsets)};

    std::vector<RawRecordBatch::RecordExtent> extents;
    for (std::size_t i{0}; i < std::size(offsets); ++i) {
        const std::uint32_t dataOffset{offsets[i]};
        const std::uint32_t blockSize{ReadU32(std::data(buffer) + dataOffset - 4)};
        extents.push_back({dataOffset, blockSize});
    }

    const RawRecordBatch batch{std::move(buffer), std::move(extents)};
    EXPECT_EQ(batch.RecordCount(), 3U);
}

TEST(RawRecordBatch, BufferSize)
{
    std::vector<std::uint32_t> offsets;
    std::vector<std::byte> buffer{BuildMultiRecordBuffer(offsets)};
    const std::size_t expectedSize{std::size(buffer)};

    std::vector<RawRecordBatch::RecordExtent> extents;
    for (std::size_t i{0}; i < std::size(offsets); ++i) {
        const std::uint32_t dataOffset{offsets[i]};
        const std::uint32_t blockSize{ReadU32(std::data(buffer) + dataOffset - 4)};
        extents.push_back({dataOffset, blockSize});
    }

    const RawRecordBatch batch{std::move(buffer), std::move(extents)};
    EXPECT_EQ(batch.BufferSize(), expectedSize);
}

}  // namespace Samoa
}  // namespace PacBio
