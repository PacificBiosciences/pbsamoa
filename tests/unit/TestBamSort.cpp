#include "TestTempDir.hpp"

#include "TestBamIo.hpp"

#include <pbsamoa/core/BamRecord.hpp>
#include <pbsamoa/core/CigarOp.hpp>
#include <pbsamoa/core/SamHeader.hpp>
#include <pbsamoa/core/Tags.hpp>
#include <pbsamoa/io/BamRawReader.hpp>
#include <pbsamoa/io/BamSort.hpp>
#include <pbsamoa/io/BamWriter.hpp>

#include <gtest/gtest.h>

#include <filesystem>
#include <format>
#include <string>
#include <thread>
#include <vector>

#include <cstdint>

#include <sys/stat.h>
#include <unistd.h>

namespace PacBio {
namespace Samoa {
namespace {

SamHeader MakeHeader(std::int32_t numRefs)
{
    SamHeader header;
    header.SetVersion("1.6");
    header.SetSortOrder("unknown");
    for (std::int32_t i{0}; i < numRefs; ++i) {
        header.AddReferenceSequence(ReferenceSequence{std::format("ref{}", i), 1000000});
    }
    return header;
}

BamRecord Mapped(std::string name, std::int32_t refId, std::int32_t pos, bool reverse)
{
    BamRecord record;
    record.Name(std::move(name))
        .Flag(reverse ? std::uint16_t{0x10} : std::uint16_t{0x0})
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

BamRecord Unmapped(std::string name, std::string sequence = "ACGT")
{
    const std::size_t sequenceLength{std::size(sequence)};
    BamRecord record;
    record.Name(std::move(name))
        .Flag(0x4)
        .RefId(-1)
        .Pos(-1)
        .MapQ(0)
        .NextRefId(-1)
        .NextPos(-1)
        .Tlen(0)
        .Sequence(std::move(sequence))
        .Qualities(std::vector<std::uint8_t>(sequenceLength, std::uint8_t{30}));
    return record;
}

// The spec coordinate key, replicated in the test so ordering assertions encode
// the documented intent (unmapped-last, forward-before-reverse).
std::uint64_t CoordKey(std::int32_t refId, std::int32_t pos, bool reverse, std::int32_t numRefs)
{
    const std::uint64_t refAdj{(refId < 0) ? static_cast<std::uint64_t>(numRefs)
                                           : static_cast<std::uint64_t>(refId)};
    const std::uint64_t posPart{static_cast<std::uint64_t>(static_cast<std::int64_t>(pos) + 1)};
    return (refAdj << 32U) | (posPart << 1U) | (reverse ? 1U : 0U);
}

}  // namespace

class BamSortTest : public ::testing::Test
{
protected:
    tests::TempDirGuard tempDir_;

    void SetUp() override { tempDir_.Reset("bam_sort"); }

    std::filesystem::path In() const { return tempDir_.File("in.bam"); }

    std::filesystem::path Out() const { return tempDir_.File("out.bam"); }
};

// --- StrNumCmp: natural / numeric-string ordering (htslib strnum_cmp parity) ---

TEST(StrNumCmpTest, EqualStrings)
{
    EXPECT_EQ(detail::StrNumCmp("read1", "read1"), 0);
    EXPECT_EQ(detail::StrNumCmp("", ""), 0);
}

TEST(StrNumCmpTest, NumericRunsCompareByValueNotLexically)
{
    // "2" < "10" numerically, even though lexically '2' > '1'.
    EXPECT_LT(detail::StrNumCmp("read2", "read10"), 0);
    EXPECT_GT(detail::StrNumCmp("read10", "read2"), 0);
    EXPECT_LT(detail::StrNumCmp("read9", "read100"), 0);
}

TEST(StrNumCmpTest, LeadingZeros)
{
    // Equal numeric value; the form with fewer characters (fewer leading zeros)
    // sorts first, matching htslib.
    EXPECT_LT(detail::StrNumCmp("a1", "a01"), 0);
    EXPECT_GT(detail::StrNumCmp("a01", "a1"), 0);
    EXPECT_LT(detail::StrNumCmp("a01", "a2"), 0);  // value 1 < 2 regardless of zeros
}

TEST(StrNumCmpTest, LexicalSegmentsAndMixed)
{
    EXPECT_LT(detail::StrNumCmp("abc", "abd"), 0);
    EXPECT_LT(detail::StrNumCmp("x9y", "x10y"), 0);  // numeric run dominates
    EXPECT_LT(detail::StrNumCmp("chr1", "chr2"), 0);
    EXPECT_GT(detail::StrNumCmp("chr11", "chr2"), 0);
    EXPECT_LT(detail::StrNumCmp("read", "read1"), 0);  // prefix sorts first
}

// --- Coordinate order: unmapped-last, forward-before-reverse, monotonic ---

TEST_F(BamSortTest, CoordinateOrderingUnmappedLastForwardBeforeReverse)
{
    std::vector<BamRecord> records;
    records.push_back(Mapped("A", 0, 100, false));
    records.push_back(Mapped("B", 0, 50, false));
    records.push_back(Mapped("C", 1, 10, false));
    records.push_back(Unmapped("D"));
    records.push_back(Mapped("E", 0, 50, true));  // same pos as B, reverse -> after B
    tests::WriteBam(In(), MakeHeader(2), records);

    const SortStats stats = SortBam(In(), Out(), SortConfig{.Order = SortOrder::COORDINATE});
    EXPECT_EQ(stats.NumRecords, 5);
    EXPECT_EQ(stats.NumRuns, 0u);  // small input -> fully in memory

    EXPECT_EQ(tests::ReadNames(Out()), (std::vector<std::string>{"B", "E", "A", "C", "D"}));
}

TEST_F(BamSortTest, StableTieBreakingPreservesInputOrder)
{
    // Ten records sharing one coordinate key must emerge in input order.
    std::vector<BamRecord> records;
    std::vector<std::string> expected;
    for (int i = 0; i < 10; ++i) {
        const std::string name{std::format("r{}", i)};
        records.push_back(Mapped(name, 0, 100, false));
        expected.push_back(name);
    }
    tests::WriteBam(In(), MakeHeader(1), records);

    SortBam(In(), Out(), SortConfig{.Order = SortOrder::COORDINATE});
    EXPECT_EQ(tests::ReadNames(Out()), expected);
}

// A FIFO is non-seekable: SortBam must read it sequentially (synchronous parsing,
// no parallel-pipeline reopen+seek, no trailing EOF-marker probe). This is what
// lets pbmm2 sort on-the-fly straight off the alignment writer's pipe. The tiny
// memory budget also proves runs spill correctly while consuming a stream.
TEST_F(BamSortTest, SortsCoordinateFromNonSeekableFifoStream)
{
    const std::filesystem::path fifo{tempDir_.File("stream.fifo")};
    ASSERT_EQ(mkfifo(fifo.c_str(), 0600), 0);

    constexpr std::int32_t NUM_RECORDS{300};
    constexpr std::int32_t NUM_REFS{3};
    std::vector<BamRecord> records;
    for (std::int32_t i{0}; i < NUM_RECORDS; ++i) {
        const std::int32_t refId{i % NUM_REFS};
        const std::int32_t pos{(i * 37) % 1000};
        records.push_back(Mapped(std::format("r{}", i), refId, pos, (i % 2) == 0));
    }
    const SamHeader header{MakeHeader(NUM_REFS)};

    // Opening the FIFO for write blocks until SortBam opens the read end, so the
    // producer must stream concurrently.
    std::thread producer{[&]() { tests::WriteBam(fifo, header, records); }};

    const SortStats stats = SortBam(
        fifo, Out(), SortConfig{.Order = SortOrder::COORDINATE, .MaxMemory = ByteLimit{1024}});
    producer.join();
    unlink(fifo.c_str());

    EXPECT_EQ(stats.NumRecords, NUM_RECORDS);
    EXPECT_GT(stats.NumRuns, 1u);  // spilled while streaming

    // The stream result must match sorting the same records from a regular file.
    tests::WriteBam(In(), header, records);
    const std::filesystem::path fileOut{tempDir_.File("file.bam")};
    SortBam(In(), fileOut, SortConfig{.Order = SortOrder::COORDINATE});
    EXPECT_EQ(tests::ReadNames(Out()), tests::ReadNames(fileOut));
}

TEST_F(BamSortTest, MultiRunMergeEqualsInMemoryAndIsMonotonic)
{
    constexpr std::int32_t NUM_RECORDS{300};
    constexpr std::int32_t NUM_REFS{3};
    std::vector<BamRecord> records;
    for (std::int32_t i{0}; i < NUM_RECORDS; ++i) {
        const std::int32_t refId{i % NUM_REFS};
        const std::int32_t pos{(i * 37) % 1000};
        records.push_back(Mapped(std::format("r{}", i), refId, pos, (i % 2) == 0));
    }
    tests::WriteBam(In(), MakeHeader(NUM_REFS), records);

    // Whole input in memory.
    const std::filesystem::path memOut{tempDir_.File("mem.bam")};
    const SortStats memStats = SortBam(In(), memOut, SortConfig{.Order = SortOrder::COORDINATE});
    EXPECT_EQ(memStats.NumRuns, 0u);

    // Tiny budget forces many spilled runs.
    const std::filesystem::path spillOut{tempDir_.File("spill.bam")};
    const SortStats spillStats = SortBam(
        In(), spillOut, SortConfig{.Order = SortOrder::COORDINATE, .MaxMemory = ByteLimit{1024}});
    EXPECT_GT(spillStats.NumRuns, 1u);

    const std::vector<std::string> memNames{tests::ReadNames(memOut)};
    const std::vector<std::string> spillNames{tests::ReadNames(spillOut)};
    EXPECT_EQ(memNames, spillNames);  // determinism: spill path == in-memory path
    EXPECT_EQ(static_cast<std::int32_t>(memNames.size()), NUM_RECORDS);

    // Output is actually coordinate-sorted (monotonic non-decreasing key).
    BamRawReader reader{spillOut};
    std::uint64_t previous{0};
    for (const auto& view : reader.Records()) {
        const std::uint64_t key{
            CoordKey(view.RefId(), view.Pos(), view.IsReverseStrand(), NUM_REFS)};
        EXPECT_GE(key, previous);
        previous = key;
    }
}

TEST_F(BamSortTest, MinimiserClustersUnmappedReadsLikeSamtoolsForwardStrand)
{
    const std::filesystem::path destinationDir{tempDir_.File("destination")};
    std::filesystem::create_directories(destinationDir);
    const std::filesystem::path memoryOutput{destinationDir / "memory.bam"};
    const std::filesystem::path spilledOutput{destinationDir / "spilled.bam"};
    tests::WriteBam(
        In(), MakeHeader(1),
        {Mapped("mapped", 0, 10, false), Unmapped("polyA", "AAAAAAAAAAAAAAAAAAAA"),
         Unmapped("likeA", "AAAAAAAAAAAAAAAAAAAC"), Unmapped("polyC", "CCCCCCCCCCCCCCCCCCCC"),
         Unmapped("mixed", "ACGTACGTACGTACGTACGT"), Unmapped("short", "ACGT")});

    const SortStats memoryStats{SortBam(
        In(), memoryOutput,
        SortConfig{.Order = SortOrder::COORDINATE, .TempDir = destinationDir, .Minimise = true})};
    const SortStats spilledStats{SortBam(In(), spilledOutput,
                                         SortConfig{.Order = SortOrder::COORDINATE,
                                                    .MaxMemory = ByteLimit{1},
                                                    .TempDir = destinationDir,
                                                    .Minimise = true})};

    EXPECT_EQ(memoryStats.NumRuns, 0U);
    EXPECT_GT(spilledStats.NumRuns, 0U);
    const std::vector<std::string> expected{"mapped", "short", "polyC", "mixed", "polyA", "likeA"};
    EXPECT_EQ(tests::ReadNames(memoryOutput), expected);
    EXPECT_EQ(tests::ReadNames(spilledOutput), expected);
    const BamRawReader reader{spilledOutput};
    EXPECT_EQ(reader.Header().SortOrder(), "coordinate");
    EXPECT_EQ(reader.Header().SubSort(), "coordinate:minhash");
}

// --- Tag order: missing-tag-first; numeric and string tag values ---

TEST_F(BamSortTest, TagSortNumericMissingFirst)
{
    const TagKey key{'x', 's'};
    auto withTag = [&](std::string name, std::int64_t value) {
        BamRecord record{Mapped(std::move(name), 0, 0, false)};  // identical coordinate
        TagMap tags;
        tags.Set(key, TagValue{value});
        record.Tags(tags);
        return record;
    };

    std::vector<BamRecord> records;
    records.push_back(withTag("v3", 3));
    records.push_back(Mapped("missing", 0, 0, false));  // no xs tag -> sorts first
    records.push_back(withTag("v1", 1));
    records.push_back(withTag("v2", 2));
    tests::WriteBam(In(), MakeHeader(1), records);

    SortBam(In(), Out(), SortConfig{.Order = SortOrder::TAG, .Tag = {'x', 's'}});
    EXPECT_EQ(tests::ReadNames(Out()), (std::vector<std::string>{"missing", "v1", "v2", "v3"}));
}

TEST_F(BamSortTest, TagSortStringValues)
{
    const TagKey key{'B', 'C'};
    auto withTag = [&](std::string name, std::string value) {
        BamRecord record{Mapped(std::move(name), 0, 0, false)};
        TagMap tags;
        tags.Set(key, TagValue{std::move(value)});
        record.Tags(tags);
        return record;
    };

    std::vector<BamRecord> records;
    records.push_back(withTag("banana", "banana"));
    records.push_back(withTag("apple", "apple"));
    records.push_back(withTag("cherry", "cherry"));
    tests::WriteBam(In(), MakeHeader(1), records);

    SortBam(In(), Out(), SortConfig{.Order = SortOrder::TAG, .Tag = {'B', 'C'}});
    EXPECT_EQ(tests::ReadNames(Out()), (std::vector<std::string>{"apple", "banana", "cherry"}));
}

// --- Zero records: header-only output with the sort order set ---

TEST_F(BamSortTest, ZeroRecordsProducesHeaderOnlyWithSortOrder)
{
    tests::WriteBam(In(), MakeHeader(1), {});

    const SortStats stats = SortBam(In(), Out(), SortConfig{.Order = SortOrder::COORDINATE});
    EXPECT_EQ(stats.NumRecords, 0);
    EXPECT_EQ(stats.NumRuns, 0u);

    BamRawReader reader{Out()};
    EXPECT_EQ(reader.Header().SortOrder(), "coordinate");
    EXPECT_FALSE(reader.ReadRecord());
}

}  // namespace Samoa
}  // namespace PacBio
