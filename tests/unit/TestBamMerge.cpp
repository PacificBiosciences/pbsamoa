#include "TestTempDir.hpp"

#include "TestBamIo.hpp"

#include <pbsamoa/core/BamRecord.hpp>
#include <pbsamoa/core/CigarOp.hpp>
#include <pbsamoa/core/SamHeader.hpp>
#include <pbsamoa/core/Tags.hpp>
#include <pbsamoa/io/BamMerge.hpp>
#include <pbsamoa/io/BamRawReader.hpp>
#include <pbsamoa/io/BamWriter.hpp>

#include <gtest/gtest.h>

#include <filesystem>
#include <format>
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

    BamRawReader reader{Out()};
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

}  // namespace Samoa
}  // namespace PacBio
