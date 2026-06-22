#include "TestTempDir.hpp"

#include "TestBamIo.hpp"

#include <pbsamoa/core/BamRecord.hpp>
#include <pbsamoa/core/CigarOp.hpp>
#include <pbsamoa/core/RawRecord.hpp>
#include <pbsamoa/core/SamHeader.hpp>
#include <pbsamoa/io/BamRawReader.hpp>
#include <pbsamoa/io/MergeReadAhead.hpp>

#include <gtest/gtest.h>

#include <filesystem>
#include <format>
#include <optional>
#include <span>
#include <string>
#include <vector>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {
namespace {

SamHeader MakeHeader(std::int32_t numRefs, std::string_view sortOrder)
{
    SamHeader header;
    header.SetVersion("1.6");
    header.SetSortOrder(std::string{sortOrder});
    for (std::int32_t i{0}; i < numRefs; ++i) {
        header.AddReferenceSequence(ReferenceSequence{std::format("ref{}", i), 1000000});
    }
    return header;
}

BamRecord Mapped(std::string name, std::int32_t refId, std::int32_t pos)
{
    BamRecord record;
    record.Name(std::move(name))
        .Flag(0)
        .RefId(refId)
        .Pos(pos)
        .MapQ(30)
        .Cigar({CigarOp{CigarOpType::M, 4}})
        .NextRefId(-1)
        .NextPos(-1)
        .Tlen(0)
        .Sequence("ACGT")
        .Qualities({30, 30, 30, 30});
    return record;
}

/// Build a coordinate-ascending BAM with \p count records. With count large
/// enough the BGZF output spans many blocks (and records straddle block
/// boundaries), exercising the read-ahead framing path.
std::vector<std::string> WriteManyRecords(const std::filesystem::path& path, std::int32_t count)
{
    std::vector<BamRecord> records;
    std::vector<std::string> names;
    records.reserve(static_cast<std::size_t>(count));
    names.reserve(static_cast<std::size_t>(count));
    for (std::int32_t i{0}; i < count; ++i) {
        std::string name{std::format("r{:06}", i)};
        names.push_back(name);
        records.push_back(Mapped(std::move(name), 0, i));
    }
    tests::WriteBam(path, MakeHeader(1, "coordinate"), records);
    return names;
}

std::vector<std::string> DrainSource(MergeReadAhead& ra, std::size_t sourceIndex)
{
    std::vector<std::string> names;
    while (const std::optional<std::span<const std::byte>> bytes{ra.Next(sourceIndex)}) {
        names.emplace_back(RawRecord{*bytes}.Name());
    }
    return names;
}

class MergeReadAheadTest : public ::testing::Test
{
protected:
    tests::TempDirGuard tempDir_;

    void SetUp() override { tempDir_.Reset("merge_readahead"); }

    std::filesystem::path File(std::string_view name) const { return tempDir_.File(name); }
};

TEST_F(MergeReadAheadTest, YieldsAllRecordsInFileOrder)
{
    const std::filesystem::path input{File("in.bam")};
    const std::vector<std::string> expected{WriteManyRecords(input, 2000)};

    const std::vector<std::filesystem::path> inputs{input};
    MergeReadAhead ra{inputs, /*decodeWorkers=*/4, ByteLimit{std::size_t{64} * 1024 * 1024}};

    ASSERT_EQ(ra.NumSources(), 1u);
    EXPECT_EQ(DrainSource(ra, 0), expected);
}

TEST_F(MergeReadAheadTest, SingleDecodeWorkerYieldsAllRecords)
{
    const std::filesystem::path input{File("in.bam")};
    const std::vector<std::string> expected{WriteManyRecords(input, 1500)};

    const std::vector<std::filesystem::path> inputs{input};
    MergeReadAhead ra{inputs, /*decodeWorkers=*/1, ByteLimit{std::size_t{64} * 1024 * 1024}};
    EXPECT_EQ(DrainSource(ra, 0), expected);
}

TEST_F(MergeReadAheadTest, MultipleSourcesYieldIndependentlyInFileOrder)
{
    const std::filesystem::path a{File("a.bam")};
    const std::filesystem::path b{File("b.bam")};
    const std::filesystem::path c{File("c.bam")};
    const std::vector<std::string> expectedA{WriteManyRecords(a, 800)};
    const std::vector<std::string> expectedB{WriteManyRecords(b, 1200)};
    const std::vector<std::string> expectedC{WriteManyRecords(c, 50)};

    const std::vector<std::filesystem::path> inputs{a, b, c};
    MergeReadAhead ra{inputs, /*decodeWorkers=*/4, ByteLimit{std::size_t{16} * 1024 * 1024}};

    ASSERT_EQ(ra.NumSources(), 3u);
    // Interleave the pulls across sources to mimic the heap consumer's access.
    EXPECT_EQ(DrainSource(ra, 2), expectedC);
    EXPECT_EQ(DrainSource(ra, 0), expectedA);
    EXPECT_EQ(DrainSource(ra, 1), expectedB);
}

TEST_F(MergeReadAheadTest, TinyBudgetStillCompletesWithoutHanging)
{
    // A budget far smaller than one record must not deadlock: each source always
    // keeps its head available, so the merge can still drain every record in order.
    const std::filesystem::path a{File("a.bam")};
    const std::filesystem::path b{File("b.bam")};
    const std::vector<std::string> expectedA{WriteManyRecords(a, 1000)};
    const std::vector<std::string> expectedB{WriteManyRecords(b, 1000)};

    const std::vector<std::filesystem::path> inputs{a, b};
    MergeReadAhead ra{inputs, /*decodeWorkers=*/4, ByteLimit{std::size_t{1}}};

    EXPECT_EQ(DrainSource(ra, 0), expectedA);
    EXPECT_EQ(DrainSource(ra, 1), expectedB);
}

TEST_F(MergeReadAheadTest, CustomBatchBytesDeliversEveryRecordInOrder)
{
    // A batch target of one byte forces a flush per record: the handoff granularity
    // must not affect record framing or per-source order.
    const std::filesystem::path a{File("a.bam")};
    const std::filesystem::path b{File("b.bam")};
    const std::vector<std::string> expectedA{WriteManyRecords(a, 700)};
    const std::vector<std::string> expectedB{WriteManyRecords(b, 300)};

    const std::vector<std::filesystem::path> inputs{a, b};
    MergeReadAhead ra{inputs, /*decodeWorkers=*/4, ByteLimit{std::size_t{64} * 1024 * 1024},
                      /*batchBytes=*/1};

    EXPECT_EQ(DrainSource(ra, 0), expectedA);
    EXPECT_EQ(DrainSource(ra, 1), expectedB);
}

TEST_F(MergeReadAheadTest, MissingInputThrowsOnConstruction)
{
    const std::vector<std::filesystem::path> inputs{File("does_not_exist.bam")};
    EXPECT_THROW((MergeReadAhead{inputs, 4, ByteLimit{std::size_t{1} << 20}}), std::runtime_error);
}

TEST_F(MergeReadAheadTest, PropagatesErrorOnTruncatedInput)
{
    const std::filesystem::path input{File("in.bam")};
    WriteManyRecords(input, 2000);

    // Lop bytes off the end so the trailing block (and EOF marker) are incomplete.
    const auto size{std::filesystem::file_size(input)};
    ASSERT_GT(size, 40u);
    std::filesystem::resize_file(input, size - 20);

    const std::vector<std::filesystem::path> inputs{input};
    MergeReadAhead ra{inputs, /*decodeWorkers=*/4, ByteLimit{std::size_t{64} * 1024 * 1024}};

    // Valid leading records may be delivered first; draining must ultimately throw.
    EXPECT_THROW(DrainSource(ra, 0), std::runtime_error);
}

}  // namespace
}  // namespace Samoa
}  // namespace PacBio
