#include "TestTempDir.hpp"

#include "TestBamIo.hpp"

#include <pbsamoa/core/BamRecord.hpp>
#include <pbsamoa/core/CigarOp.hpp>
#include <pbsamoa/core/SamHeader.hpp>
#include <pbsamoa/core/Tags.hpp>
#include <pbsamoa/index/BaiIndex.hpp>
#include <pbsamoa/io/BamMerge.hpp>
#include <pbsamoa/io/BamRawReader.hpp>
#include <pbsamoa/io/BamWriter.hpp>

#include <gtest/gtest.h>

#include <algorithm>
#include <filesystem>
#include <format>
#include <fstream>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

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

}  // namespace

class BamMergeTest : public ::testing::Test
{
protected:
    tests::TempDirGuard tempDir_;

    void SetUp() override { tempDir_.Reset("bam_merge"); }

    std::filesystem::path A() const { return tempDir_.File("a.bam"); }

    std::filesystem::path B() const { return tempDir_.File("b.bam"); }

    std::filesystem::path C() const { return tempDir_.File("c.bam"); }

    std::filesystem::path Out() const { return tempDir_.File("out.bam"); }
};

TEST_F(BamMergeTest, CoordinateMergeInterleavesSortedInputs)
{
    // Each input is individually coordinate-sorted; the merge interleaves them.
    tests::WriteBam(A(), MakeHeader(1, "coordinate"), {Mapped("a10", 0, 10), Mapped("a50", 0, 50)});
    tests::WriteBam(B(), MakeHeader(1, "coordinate"), {Mapped("b20", 0, 20), Mapped("b60", 0, 60)});

    const MergeStats stats =
        MergeBam({A(), B()}, Out(), MergeConfig{.Order = SortOrder::COORDINATE});
    EXPECT_EQ(stats.NumRecords, 4);
    EXPECT_EQ(stats.NumInputs, 2u);

    EXPECT_EQ(tests::ReadNames(Out()), (std::vector<std::string>{"a10", "b20", "a50", "b60"}));

    const BamRawReader reader{Out()};
    EXPECT_EQ(reader.Header().SortOrder(), "coordinate");
}

TEST_F(BamMergeTest, CoordinateMergeAcrossThreeInputsAndReferences)
{
    tests::WriteBam(A(), MakeHeader(2, "coordinate"), {Mapped("a", 0, 5), Mapped("a2", 1, 5)});
    tests::WriteBam(B(), MakeHeader(2, "coordinate"), {Mapped("b", 0, 7), Mapped("b2", 1, 1)});
    tests::WriteBam(C(), MakeHeader(2, "coordinate"), {Mapped("c", 0, 1)});

    MergeBam({A(), B(), C()}, Out(), MergeConfig{.Order = SortOrder::COORDINATE});

    // ref0: pos1(c), pos5(a), pos7(b); then ref1: pos1(b2), pos5(a2)
    EXPECT_EQ(tests::ReadNames(Out()), (std::vector<std::string>{"c", "a", "b", "b2", "a2"}));
}

TEST_F(BamMergeTest, QueryNameMergeUsesNaturalOrder)
{
    tests::WriteBam(A(), MakeHeader(1, "queryname"), {Mapped("r1", 0, 0), Mapped("r3", 0, 0)});
    tests::WriteBam(B(), MakeHeader(1, "queryname"), {Mapped("r2", 0, 0), Mapped("r10", 0, 0)});

    MergeBam({A(), B()}, Out(), MergeConfig{.Order = SortOrder::QUERY_NAME});

    // Natural order: r1 < r2 < r3 < r10 (numeric, not lexical).
    EXPECT_EQ(tests::ReadNames(Out()), (std::vector<std::string>{"r1", "r2", "r3", "r10"}));
}

TEST_F(BamMergeTest, TagMergeMissingFirst)
{
    const TagKey key{'x', 's'};
    auto withTag = [&](std::string name, std::int64_t value) {
        BamRecord record{Mapped(std::move(name), 0, 0)};
        TagMap tags;
        tags.Set(key, TagValue{value});
        record.Tags(tags);
        return record;
    };

    // Tag order inputs carry SO:unknown; each is independently tag-sorted.
    tests::WriteBam(A(), MakeHeader(1, "unknown"), {Mapped("missing", 0, 0), withTag("a1", 1)});
    tests::WriteBam(B(), MakeHeader(1, "unknown"), {withTag("b2", 2), withTag("b3", 3)});

    MergeBam({A(), B()}, Out(), MergeConfig{.Order = SortOrder::TAG, .Tag = {'x', 's'}});

    EXPECT_EQ(tests::ReadNames(Out()), (std::vector<std::string>{"missing", "a1", "b2", "b3"}));
}

TEST_F(BamMergeTest, EqualKeyTiesBreakOnInputOrder)
{
    // Both records share one coordinate key; the first input must win the tie.
    tests::WriteBam(A(), MakeHeader(1, "coordinate"), {Mapped("from_a", 0, 100)});
    tests::WriteBam(B(), MakeHeader(1, "coordinate"), {Mapped("from_b", 0, 100)});

    MergeBam({A(), B()}, Out(), MergeConfig{.Order = SortOrder::COORDINATE});
    EXPECT_EQ(tests::ReadNames(Out()), (std::vector<std::string>{"from_a", "from_b"}));
}

TEST_F(BamMergeTest, SingleInputMergeCopies)
{
    tests::WriteBam(A(), MakeHeader(1, "coordinate"), {Mapped("x", 0, 1), Mapped("y", 0, 2)});

    const MergeStats stats = MergeBam({A()}, Out(), MergeConfig{.Order = SortOrder::COORDINATE});
    EXPECT_EQ(stats.NumRecords, 2);
    EXPECT_EQ(stats.NumInputs, 1u);
    EXPECT_EQ(tests::ReadNames(Out()), (std::vector<std::string>{"x", "y"}));
}

// --- parallel read-ahead merge ---

TEST_F(BamMergeTest, ParallelMergeIsDeterministicAcrossThreadsAndMemory)
{
    // Four coordinate-sorted inputs whose positions interleave, plus enough records
    // to span many BGZF blocks and exercise the read-ahead producers.
    constexpr std::int32_t PER_INPUT{600};
    std::vector<std::filesystem::path> inputs{A(), B(), C(), tempDir_.File("d.bam")};
    for (std::size_t s{0}; s < std::size(inputs); ++s) {
        std::vector<BamRecord> records;
        records.reserve(PER_INPUT);
        for (std::int32_t i{0}; i < PER_INPUT; ++i) {
            const std::int32_t pos{(i * 4) + static_cast<std::int32_t>(s)};
            records.push_back(Mapped(std::format("s{}_r{:04}", s, i), 0, pos));
        }
        tests::WriteBam(inputs[s], MakeHeader(1, "coordinate"), records);
    }

    const auto run{[&](std::size_t threads, ByteLimit mem) {
        MergeBam(
            inputs, Out(),
            MergeConfig{
                .Order = SortOrder::COORDINATE, .NumThreads = threads, .ReadAheadMemory = mem});
        return tests::ReadNames(Out());
    }};

    const std::vector<std::string> reference{run(1, ByteLimit{std::size_t{64} * 1024 * 1024})};
    EXPECT_EQ(std::size(reference), static_cast<std::size_t>(PER_INPUT) * std::size(inputs));
    // Many threads and a sub-record budget must reproduce the single-thread output.
    EXPECT_EQ(run(8, ByteLimit{std::size_t{64} * 1024 * 1024}), reference);
    EXPECT_EQ(run(8, ByteLimit{std::size_t{1}}), reference);
}

TEST_F(BamMergeTest, GranularThreadAndMemoryKnobsReproduceDefaultOutput)
{
    // The decode/compress/batch/writer-queue knobs change only resource usage, not
    // the merged output: every combination must reproduce the default-config result.
    constexpr std::int32_t PER_INPUT{400};
    const std::vector<std::filesystem::path> inputs{A(), B(), C()};
    for (std::size_t s{0}; s < std::size(inputs); ++s) {
        std::vector<BamRecord> records;
        records.reserve(PER_INPUT);
        for (std::int32_t i{0}; i < PER_INPUT; ++i) {
            records.push_back(
                Mapped(std::format("s{}_r{:04}", s, i), 0, (i * 3) + static_cast<std::int32_t>(s)));
        }
        tests::WriteBam(inputs[s], MakeHeader(1, "coordinate"), records);
    }

    const auto run{[&](const MergeConfig& cfg) {
        MergeBam(inputs, Out(), cfg);
        return tests::ReadNames(Out());
    }};

    const std::vector<std::string> reference{run(MergeConfig{.Order = SortOrder::COORDINATE})};
    EXPECT_EQ(std::size(reference), static_cast<std::size_t>(PER_INPUT) * std::size(inputs));

    // Asymmetric pools, a sub-record batch, and a one-deep writer queue.
    EXPECT_EQ(run(MergeConfig{.Order = SortOrder::COORDINATE,
                              .DecodeThreads = 4,
                              .CompressThreads = 1,
                              .BatchBytes = 128,
                              .WriterQueueCapacity = 1}),
              reference);
    // Decode/compress inherit NumThreads at 0; a batch larger than the whole stream.
    EXPECT_EQ(run(MergeConfig{.Order = SortOrder::COORDINATE,
                              .NumThreads = 3,
                              .BatchBytes = std::size_t{8} * 1024 * 1024,
                              .WriterQueueCapacity = 512}),
              reference);
}

TEST_F(BamMergeTest, SkewedInputsMergeInCoordinateOrder)
{
    // A holds the low coordinates, B the high ones: the merge drains A entirely
    // before B contributes, exercising backpressure on the idle source.
    std::vector<BamRecord> lo;
    std::vector<BamRecord> hi;
    std::vector<std::string> expected;
    for (std::int32_t i{0}; i < 400; ++i) {
        lo.push_back(Mapped(std::format("lo{:04}", i), 0, i));
        expected.push_back(std::format("lo{:04}", i));
    }
    for (std::int32_t i{0}; i < 400; ++i) {
        hi.push_back(Mapped(std::format("hi{:04}", i), 0, 100000 + i));
    }
    for (std::int32_t i{0}; i < 400; ++i) {
        expected.push_back(std::format("hi{:04}", i));
    }
    tests::WriteBam(A(), MakeHeader(1, "coordinate"), lo);
    tests::WriteBam(B(), MakeHeader(1, "coordinate"), hi);

    MergeBam({A(), B()}, Out(),
             MergeConfig{.Order = SortOrder::COORDINATE,
                         .NumThreads = 8,
                         .ReadAheadMemory = ByteLimit{std::size_t{4} * 1024}});
    EXPECT_EQ(tests::ReadNames(Out()), expected);
}

TEST_F(BamMergeTest, IncompatibleReferencesThrows)
{
    tests::WriteBam(A(), MakeHeader(1, "coordinate"), {Mapped("a", 0, 1)});

    // Different reference name/length than A.
    SamHeader headerB;
    headerB.SetVersion("1.6");
    headerB.SetSortOrder("coordinate");
    headerB.AddReferenceSequence(ReferenceSequence{"different", 500});
    tests::WriteBam(B(), headerB, {Mapped("b", 0, 1)});

    EXPECT_THROW(MergeBam({A(), B()}, Out(), MergeConfig{.Order = SortOrder::COORDINATE}),
                 std::runtime_error);
}

TEST_F(BamMergeTest, SortOrderMismatchThrows)
{
    tests::WriteBam(A(), MakeHeader(1, "coordinate"), {Mapped("a", 0, 1)});
    tests::WriteBam(B(), MakeHeader(1, "coordinate"), {Mapped("b", 0, 2)});

    // Inputs are coordinate-sorted but a queryname merge is requested.
    EXPECT_THROW(MergeBam({A(), B()}, Out(), MergeConfig{.Order = SortOrder::QUERY_NAME}),
                 std::runtime_error);
}

TEST_F(BamMergeTest, EmptyInputListThrows)
{
    EXPECT_THROW(MergeBam({}, Out(), MergeConfig{}), std::runtime_error);
}

// --- auto-detection ---

TEST_F(BamMergeTest, AutoDetectsCoordinateMerge)
{
    tests::WriteBam(A(), MakeHeader(1, "coordinate"), {Mapped("a10", 0, 10), Mapped("a50", 0, 50)});
    tests::WriteBam(B(), MakeHeader(1, "coordinate"), {Mapped("b20", 0, 20), Mapped("b60", 0, 60)});

    // No explicit order: detected from @HD SO:coordinate and interleaved.
    const MergeStats stats = MergeBam({A(), B()}, Out(), MergeConfig{});
    EXPECT_EQ(stats.NumRecords, 4);
    EXPECT_EQ(tests::ReadNames(Out()), (std::vector<std::string>{"a10", "b20", "a50", "b60"}));

    const BamRawReader reader{Out()};
    EXPECT_EQ(reader.Header().SortOrder(), "coordinate");
}

TEST_F(BamMergeTest, AutoDetectsQueryNameMerge)
{
    tests::WriteBam(A(), MakeHeader(1, "queryname"), {Mapped("r1", 0, 0), Mapped("r3", 0, 0)});
    tests::WriteBam(B(), MakeHeader(1, "queryname"), {Mapped("r2", 0, 0), Mapped("r10", 0, 0)});

    MergeBam({A(), B()}, Out(), MergeConfig{});
    EXPECT_EQ(tests::ReadNames(Out()), (std::vector<std::string>{"r1", "r2", "r3", "r10"}));
}

TEST_F(BamMergeTest, MixedSortOrdersThrows)
{
    tests::WriteBam(A(), MakeHeader(1, "coordinate"), {Mapped("a", 0, 1)});
    tests::WriteBam(B(), MakeHeader(1, "queryname"), {Mapped("b", 0, 2)});

    // Cannot auto-detect a single order across coordinate + queryname inputs.
    EXPECT_THROW(MergeBam({A(), B()}, Out(), MergeConfig{}), std::runtime_error);
}

TEST_F(BamMergeTest, MixedSortedAndUnsortedThrows)
{
    tests::WriteBam(A(), MakeHeader(1, "coordinate"), {Mapped("a", 0, 1)});
    tests::WriteBam(B(), MakeHeader(1, "unsorted"), {Mapped("b", 0, 2)});

    EXPECT_THROW(MergeBam({A(), B()}, Out(), MergeConfig{}), std::runtime_error);
}

// --- concat (block passthrough) ---

TEST_F(BamMergeTest, AutoDetectsConcatForUnsortedInputs)
{
    tests::WriteBam(A(), MakeHeader(1, "unsorted"), {Mapped("a1", 0, 50), Mapped("a2", 0, 10)});
    tests::WriteBam(B(), MakeHeader(1, "unsorted"), {Mapped("b1", 0, 30)});

    const MergeStats stats = MergeBam({A(), B()}, Out(), MergeConfig{});
    // Concat preserves per-file order and appends files: no interleaving.
    EXPECT_EQ(tests::ReadNames(Out()), (std::vector<std::string>{"a1", "a2", "b1"}));
    EXPECT_EQ(stats.NumInputs, 2u);
    EXPECT_EQ(stats.NumRecords, -1);  // not counted in passthrough

    const BamRawReader reader{Out()};
    EXPECT_EQ(reader.Header().SortOrder(), "unsorted");
}

TEST_F(BamMergeTest, ConcatFlagIgnoresSortOrder)
{
    // Inputs are coordinate-sorted, but --concat forces byte concatenation: the
    // output is the inputs appended in order, not a coordinate interleave.
    tests::WriteBam(A(), MakeHeader(1, "coordinate"), {Mapped("a10", 0, 10), Mapped("a50", 0, 50)});
    tests::WriteBam(B(), MakeHeader(1, "coordinate"), {Mapped("b20", 0, 20), Mapped("b60", 0, 60)});

    MergeBam({A(), B()}, Out(), MergeConfig{.Concat = true});
    EXPECT_EQ(tests::ReadNames(Out()), (std::vector<std::string>{"a10", "a50", "b20", "b60"}));

    const BamRawReader reader{Out()};
    EXPECT_EQ(reader.Header().SortOrder(), "unsorted");
}

TEST_F(BamMergeTest, ConcatUnionsHeaderComments)
{
    SamHeader headerA{MakeHeader(1, "unsorted")};
    headerA.AddComment("from-A");
    SamHeader headerB{MakeHeader(1, "unsorted")};
    headerB.AddComment("from-B");
    tests::WriteBam(A(), headerA, {Mapped("a", 0, 1)});
    tests::WriteBam(B(), headerB, {Mapped("b", 0, 1)});

    MergeBam({A(), B()}, Out(), MergeConfig{.Concat = true});

    const BamRawReader reader{Out()};
    const std::span<const std::string> comments{reader.Header().Comments()};
    EXPECT_NE(std::ranges::find(comments, "from-A"), std::end(comments));
    EXPECT_NE(std::ranges::find(comments, "from-B"), std::end(comments));
}

TEST_F(BamMergeTest, ConcatIncompatibleReferencesThrows)
{
    tests::WriteBam(A(), MakeHeader(1, "unsorted"), {Mapped("a", 0, 1)});

    SamHeader headerB;
    headerB.SetVersion("1.6");
    headerB.SetSortOrder("unsorted");
    headerB.AddReferenceSequence(ReferenceSequence{"different", 500});
    tests::WriteBam(B(), headerB, {Mapped("b", 0, 1)});

    EXPECT_THROW(MergeBam({A(), B()}, Out(), MergeConfig{.Concat = true}), std::runtime_error);
}

// --- parallelism determinism ---

TEST_F(BamMergeTest, CoordinateMergeDeterministicAcrossThreadCounts)
{
    tests::WriteBam(A(), MakeHeader(2, "coordinate"),
                    {Mapped("a", 0, 5), Mapped("a2", 1, 5), Mapped("a3", 0, 90)});
    tests::WriteBam(B(), MakeHeader(2, "coordinate"),
                    {Mapped("b", 0, 7), Mapped("b2", 1, 1), Mapped("b3", 0, 5)});
    tests::WriteBam(C(), MakeHeader(2, "coordinate"), {Mapped("c", 0, 1)});

    const std::filesystem::path out1{tempDir_.File("t1.bam")};
    const std::filesystem::path out8{tempDir_.File("t8.bam")};
    MergeBam({A(), B(), C()}, out1, MergeConfig{.Order = SortOrder::COORDINATE, .NumThreads = 1});
    MergeBam({A(), B(), C()}, out8, MergeConfig{.Order = SortOrder::COORDINATE, .NumThreads = 8});

    // Record order (and equal-key tie-breaking) must not depend on the thread count.
    EXPECT_EQ(tests::ReadNames(out1), tests::ReadNames(out8));
}

// --- disjoint-chain verbatim passthrough ---
//
// A coordinate merge whose inputs occupy strictly disjoint coordinate ranges must
// emit via verbatim BGZF block copy (no recompression) rather than the heap merge.
// MergeStats::Passthrough makes the chosen path observable so these tests pin the
// optimization, not merely the (identical) record output.

TEST_F(BamMergeTest, DisjointCoordinateInputsUsePassthrough)
{
    tests::WriteBam(A(), MakeHeader(1, "coordinate"), {Mapped("a1", 0, 10), Mapped("a2", 0, 20)});
    tests::WriteBam(B(), MakeHeader(1, "coordinate"), {Mapped("b1", 0, 30), Mapped("b2", 0, 40)});

    const MergeStats stats =
        MergeBam({A(), B()}, Out(), MergeConfig{.Order = SortOrder::COORDINATE});

    EXPECT_TRUE(stats.Passthrough);
    EXPECT_EQ(stats.NumRecords, 4);  // counted exactly even on the passthrough path
    EXPECT_EQ(stats.NumInputs, 2u);
    EXPECT_EQ(tests::ReadNames(Out()), (std::vector<std::string>{"a1", "a2", "b1", "b2"}));

    const BamRawReader reader{Out()};
    EXPECT_EQ(reader.Header().SortOrder(), "coordinate");
}

TEST_F(BamMergeTest, DisjointInputsPassedOutOfMinOrderAreReordered)
{
    // High-coordinate file given first; the merge must reorder by minimum coordinate.
    tests::WriteBam(A(), MakeHeader(1, "coordinate"),
                    {Mapped("hi1", 0, 300), Mapped("hi2", 0, 400)});
    tests::WriteBam(B(), MakeHeader(1, "coordinate"), {Mapped("lo1", 0, 10), Mapped("lo2", 0, 20)});

    const MergeStats stats =
        MergeBam({A(), B()}, Out(), MergeConfig{.Order = SortOrder::COORDINATE});

    EXPECT_TRUE(stats.Passthrough);
    EXPECT_EQ(tests::ReadNames(Out()), (std::vector<std::string>{"lo1", "lo2", "hi1", "hi2"}));
}

TEST_F(BamMergeTest, DisjointChainAcrossReferencesUsesPassthrough)
{
    // Each file lives on its own reference: ref0 keys all precede ref1 keys.
    tests::WriteBam(A(), MakeHeader(2, "coordinate"), {Mapped("r1a", 1, 5), Mapped("r1b", 1, 9)});
    tests::WriteBam(B(), MakeHeader(2, "coordinate"), {Mapped("r0a", 0, 5), Mapped("r0b", 0, 9)});

    const MergeStats stats =
        MergeBam({A(), B()}, Out(), MergeConfig{.Order = SortOrder::COORDINATE});

    EXPECT_TRUE(stats.Passthrough);
    EXPECT_EQ(tests::ReadNames(Out()), (std::vector<std::string>{"r0a", "r0b", "r1a", "r1b"}));
}

TEST_F(BamMergeTest, ThreeDisjointInputsChainReorderedAndCounted)
{
    tests::WriteBam(A(), MakeHeader(1, "coordinate"), {Mapped("mid", 0, 100)});
    tests::WriteBam(B(), MakeHeader(1, "coordinate"),
                    {Mapped("hi1", 0, 200), Mapped("hi2", 0, 250)});
    tests::WriteBam(C(), MakeHeader(1, "coordinate"), {Mapped("lo1", 0, 1), Mapped("lo2", 0, 2)});

    const MergeStats stats =
        MergeBam({A(), B(), C()}, Out(), MergeConfig{.Order = SortOrder::COORDINATE});

    EXPECT_TRUE(stats.Passthrough);
    EXPECT_EQ(stats.NumRecords, 5);
    EXPECT_EQ(tests::ReadNames(Out()),
              (std::vector<std::string>{"lo1", "lo2", "mid", "hi1", "hi2"}));
}

TEST_F(BamMergeTest, OverlappingInputsFallBackToHeap)
{
    // Ranges overlap (A:[10,50], B:[20,60]); passthrough is unsafe, heap must run.
    tests::WriteBam(A(), MakeHeader(1, "coordinate"), {Mapped("a10", 0, 10), Mapped("a50", 0, 50)});
    tests::WriteBam(B(), MakeHeader(1, "coordinate"), {Mapped("b20", 0, 20), Mapped("b60", 0, 60)});

    const MergeStats stats =
        MergeBam({A(), B()}, Out(), MergeConfig{.Order = SortOrder::COORDINATE});

    EXPECT_FALSE(stats.Passthrough);
    EXPECT_EQ(stats.NumRecords, 4);
    EXPECT_EQ(tests::ReadNames(Out()), (std::vector<std::string>{"a10", "b20", "a50", "b60"}));
}

TEST_F(BamMergeTest, EqualBoundaryKeyFallsBackToHeap)
{
    // max(A) == min(B) at (ref0,pos100): not strictly disjoint, so the cross-file
    // equal-key tie-break must decide — only the heap path honors it.
    tests::WriteBam(A(), MakeHeader(1, "coordinate"),
                    {Mapped("a10", 0, 10), Mapped("a100", 0, 100)});
    tests::WriteBam(B(), MakeHeader(1, "coordinate"),
                    {Mapped("b100", 0, 100), Mapped("b200", 0, 200)});

    const MergeStats stats =
        MergeBam({A(), B()}, Out(), MergeConfig{.Order = SortOrder::COORDINATE});

    EXPECT_FALSE(stats.Passthrough);
    EXPECT_EQ(tests::ReadNames(Out()), (std::vector<std::string>{"a10", "a100", "b100", "b200"}));
}

TEST_F(BamMergeTest, EmptyInputAmongDisjointStillPassesThrough)
{
    tests::WriteBam(A(), MakeHeader(1, "coordinate"), {Mapped("lo", 0, 10)});
    tests::WriteBam(B(), MakeHeader(1, "coordinate"), {});  // header only, no records
    tests::WriteBam(C(), MakeHeader(1, "coordinate"), {Mapped("hi", 0, 90)});

    const MergeStats stats =
        MergeBam({A(), B(), C()}, Out(), MergeConfig{.Order = SortOrder::COORDINATE});

    EXPECT_TRUE(stats.Passthrough);
    EXPECT_EQ(stats.NumRecords, 2);
    EXPECT_EQ(tests::ReadNames(Out()), (std::vector<std::string>{"lo", "hi"}));
}

TEST_F(BamMergeTest, QueryNameDisjointNeverPassesThrough)
{
    // Even with non-overlapping names, only coordinate merges may pass through.
    tests::WriteBam(A(), MakeHeader(1, "queryname"), {Mapped("r1", 0, 0), Mapped("r2", 0, 0)});
    tests::WriteBam(B(), MakeHeader(1, "queryname"), {Mapped("r3", 0, 0), Mapped("r4", 0, 0)});

    const MergeStats stats =
        MergeBam({A(), B()}, Out(), MergeConfig{.Order = SortOrder::QUERY_NAME});

    EXPECT_FALSE(stats.Passthrough);
    EXPECT_EQ(tests::ReadNames(Out()), (std::vector<std::string>{"r1", "r2", "r3", "r4"}));
}

namespace {

std::vector<char> ReadAllBytes(const std::filesystem::path& path)
{
    std::ifstream in{path, std::ios::binary | std::ios::ate};
    const std::streamsize size{in.tellg()};
    in.seekg(0);
    std::vector<char> bytes(static_cast<std::size_t>(size));
    in.read(std::data(bytes), size);
    return bytes;
}

}  // namespace

TEST_F(BamMergeTest, CoordinateMergeBaiMatchesBuildOnOutput)
{
    // Many records over two references so the merged output spans several BGZF blocks
    // (records cross block boundaries). The inputs' coordinate ranges overlap, forcing
    // the heap-merge path (not disjoint passthrough), which is where the index streams.
    std::vector<BamRecord> aRecs;
    std::vector<BamRecord> bRecs;
    for (std::int32_t ref{0}; ref < 2; ++ref) {
        for (std::int32_t i{0}; i < 1500; ++i) {
            aRecs.push_back(Mapped(std::format("a-{}-{}", ref, i), ref, 2 * i));
            bRecs.push_back(Mapped(std::format("b-{}-{}", ref, i), ref, (2 * i) + 1));
        }
    }
    tests::WriteBam(A(), MakeHeader(2, "coordinate"), aRecs);
    tests::WriteBam(B(), MakeHeader(2, "coordinate"), bRecs);

    const std::filesystem::path streamedBai{tempDir_.File("out.bam.bai")};
    const MergeStats stats{MergeBam(
        {A(), B()}, Out(), MergeConfig{.Order = SortOrder::COORDINATE, .BaiOutput = streamedBai})};
    ASSERT_FALSE(stats.Passthrough) << "expected heap merge, not passthrough";
    ASSERT_TRUE(std::filesystem::exists(streamedBai));

    // Oracle: scan the finished output with the file-based builder.
    const std::filesystem::path scannedBai{tempDir_.File("scanned.bai")};
    BaiIndex::Build(Out()).ToFile(scannedBai);

    EXPECT_EQ(ReadAllBytes(streamedBai), ReadAllBytes(scannedBai));
}

TEST_F(BamMergeTest, BaiOutputNotPublishedWhenBamCommitFails)
{
    tests::WriteBam(A(), MakeHeader(1, "coordinate"), {Mapped("a10", 0, 10), Mapped("a30", 0, 30)});
    tests::WriteBam(B(), MakeHeader(1, "coordinate"), {Mapped("b20", 0, 20), Mapped("b40", 0, 40)});

    ASSERT_TRUE(std::filesystem::create_directory(Out()));
    {
        std::ofstream blocker{Out() / "blocker"};
        blocker << "block";
    }

    const std::filesystem::path bai{tempDir_.File("out.bam.bai")};
    EXPECT_THROW(
        MergeBam({A(), B()}, Out(), MergeConfig{.Order = SortOrder::COORDINATE, .BaiOutput = bai}),
        std::runtime_error);
    EXPECT_FALSE(std::filesystem::exists(bai));
}

TEST_F(BamMergeTest, BaiOutputRejectedForConcat)
{
    // Concat output is SO:unsorted, which a BAI cannot index.
    tests::WriteBam(A(), MakeHeader(1, "coordinate"), {Mapped("a", 0, 10)});
    tests::WriteBam(B(), MakeHeader(1, "coordinate"), {Mapped("b", 0, 20)});

    const std::filesystem::path bai{tempDir_.File("out.bam.bai")};
    EXPECT_THROW(MergeBam({A(), B()}, Out(), MergeConfig{.Concat = true, .BaiOutput = bai}),
                 std::runtime_error);
    EXPECT_FALSE(std::filesystem::exists(bai));
}

TEST_F(BamMergeTest, BaiOutputRejectedForQuerynameOrder)
{
    tests::WriteBam(A(), MakeHeader(1, "queryname"), {Mapped("r1", 0, 0)});
    tests::WriteBam(B(), MakeHeader(1, "queryname"), {Mapped("r2", 0, 0)});

    const std::filesystem::path bai{tempDir_.File("out.bam.bai")};
    EXPECT_THROW(
        MergeBam({A(), B()}, Out(), MergeConfig{.Order = SortOrder::QUERY_NAME, .BaiOutput = bai}),
        std::runtime_error);
    EXPECT_FALSE(std::filesystem::exists(bai));
}

TEST_F(BamMergeTest, BaiOutputRejectedForDisjointPassthrough)
{
    // Disjoint coordinate ranges take the verbatim-block passthrough, which never decodes
    // records, so the on-the-fly index has nothing to observe.
    tests::WriteBam(A(), MakeHeader(1, "coordinate"), {Mapped("a1", 0, 10), Mapped("a2", 0, 20)});
    tests::WriteBam(B(), MakeHeader(1, "coordinate"), {Mapped("b1", 0, 30), Mapped("b2", 0, 40)});

    const std::filesystem::path bai{tempDir_.File("out.bam.bai")};
    EXPECT_THROW(
        MergeBam({A(), B()}, Out(), MergeConfig{.Order = SortOrder::COORDINATE, .BaiOutput = bai}),
        std::runtime_error);
    EXPECT_FALSE(std::filesystem::exists(bai));
}

}  // namespace Samoa
}  // namespace PacBio
