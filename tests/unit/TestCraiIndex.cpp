#include "TestData.hpp"

#include <pbsamoa/cram/CramCompression.hpp>
#include <pbsamoa/index/CraiIndex.hpp>

#include <gtest/gtest.h>

#include <chrono>
#include <filesystem>
#include <format>
#include <fstream>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace PacBio {
namespace Samoa {
namespace {

std::filesystem::path TempPath(std::string_view tag)
{
    const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
    return std::filesystem::temp_directory_path() / std::format("pbsamoa_{}_{}.crai", tag, stamp);
}

void WriteBinary(const std::filesystem::path& path, std::span<const std::byte> bytes)
{
    std::ofstream out{path, std::ios::binary};
    if (!out.is_open()) {
        throw std::runtime_error{"Failed to create temporary file: " + path.string()};
    }
    out.write(reinterpret_cast<const char*>(bytes.data()),
              static_cast<std::streamsize>(bytes.size()));
    if (!out.good()) {
        throw std::runtime_error{"Failed to write temporary file: " + path.string()};
    }
}

std::filesystem::path WriteGzipText(std::string_view tag, std::string_view text)
{
    const std::span<const std::byte> payload{reinterpret_cast<const std::byte*>(text.data()),
                                             std::size(text)};
    const std::vector<std::byte> compressed = CramGzipCompress(payload);
    const std::filesystem::path path = TempPath(tag);
    WriteBinary(path, compressed);
    return path;
}

class TempFileGuard
{
public:
    explicit TempFileGuard(std::filesystem::path path) : path_{std::move(path)} {}

    TempFileGuard(const TempFileGuard&) = delete;
    TempFileGuard& operator=(const TempFileGuard&) = delete;
    TempFileGuard(TempFileGuard&&) = delete;
    TempFileGuard& operator=(TempFileGuard&&) = delete;

    const std::filesystem::path& Path() const { return path_; }

    ~TempFileGuard()
    {
        std::error_code ec;
        std::filesystem::remove(path_, ec);
    }

private:
    std::filesystem::path path_;
};

}  // namespace

TEST(CraiIndex, ParsesSamtoolsIndexFixture)
{
    // Equivalent CRAI rows to the historical samtools fixture values used by this test.
    const TempFileGuard craiPath{WriteGzipText("crai_index_fixture",
                                               "0\t999901\t497\t866\t336\t8553\n"
                                               "1\t12345\t100\t9777\t123\t456\n"
                                               "2\t20000\t50\t12000\t300\t400\n"
                                               "-1\t0\t1\t14169\t346\t17346\n")};

    const CraiIndex index = CraiIndex::FromFile(craiPath.Path());
    ASSERT_EQ(std::size(index.Entries()), 4U);

    const CraiEntry& first = index.Entries().at(0);
    EXPECT_EQ(first.SequenceId, 0);
    EXPECT_EQ(first.AlignmentStart, 999901);
    EXPECT_EQ(first.AlignmentSpan, 497);
    EXPECT_EQ(first.ContainerOffset, 866);
    EXPECT_EQ(first.SliceOffset, 336);
    EXPECT_EQ(first.SliceSize, 8553);

    const CraiEntry& last = index.Entries().at(3);
    EXPECT_EQ(last.SequenceId, -1);
    EXPECT_EQ(last.AlignmentStart, 0);
    EXPECT_EQ(last.AlignmentSpan, 1);
    EXPECT_EQ(last.ContainerOffset, 14169);
    EXPECT_EQ(last.SliceOffset, 346);
    EXPECT_EQ(last.SliceSize, 17346);
}

TEST(CraiIndex, PreservesMultiReferenceDuplicateOffsets)
{
    const TempFileGuard craiPath{WriteGzipText("crai_range_fixture",
                                               "0\t100\t10\t5000\t42\t700\n"
                                               "1\t200\t20\t5000\t42\t700\n"
                                               "0\t300\t30\t6000\t50\t800\n"
                                               "1\t400\t40\t7000\t60\t900\n"
                                               "2\t500\t50\t8000\t70\t1000\n"
                                               "-1\t0\t1\t9000\t80\t1100\n")};

    const CraiIndex index = CraiIndex::FromFile(craiPath.Path());
    ASSERT_EQ(std::size(index.Entries()), 6U);

    const CraiEntry& first = index.Entries().at(0);
    const CraiEntry& second = index.Entries().at(1);
    EXPECT_EQ(first.SequenceId, 0);
    EXPECT_EQ(second.SequenceId, 1);
    EXPECT_EQ(first.ContainerOffset, second.ContainerOffset);
    EXPECT_EQ(first.SliceOffset, second.SliceOffset);
    EXPECT_EQ(first.SliceSize, second.SliceSize);
}

TEST(CraiIndex, EntriesForReferenceHandlesMappedAndUnmapped)
{
    const TempFileGuard indexPath{WriteGzipText("crai_entries_for_reference",
                                                "0\t999901\t497\t866\t336\t8553\n"
                                                "1\t12345\t100\t9777\t123\t456\n"
                                                "2\t20000\t50\t12000\t300\t400\n"
                                                "-1\t0\t1\t14169\t346\t17346\n")};
    const CraiIndex index = CraiIndex::FromFile(indexPath.Path());

    const std::vector<CraiEntry> ref1Entries = index.EntriesForReference(1);
    ASSERT_EQ(std::size(ref1Entries), 1U);
    EXPECT_EQ(ref1Entries[0].SequenceId, 1);
    EXPECT_EQ(ref1Entries[0].ContainerOffset, 9777);

    const std::vector<CraiEntry> unmappedEntries = index.EntriesForReference(-1);
    ASSERT_EQ(std::size(unmappedEntries), 1U);
    EXPECT_EQ(unmappedEntries[0].SequenceId, -1);
    EXPECT_EQ(unmappedEntries[0].SliceOffset, 346);

    EXPECT_TRUE(index.EntriesForReference(999).empty());
}

TEST(CraiIndex, ThrowsOnWrongColumnCount)
{
    const TempFileGuard input{WriteGzipText("crai_wrong_columns", "0\t1\t2\t3\t4\n")};

    try {
        (void)CraiIndex::FromFile(input.Path());
        FAIL() << "Expected parsing to throw";
    } catch (const std::runtime_error& e) {
        const std::string message = e.what();
        EXPECT_NE(message.find("CraiIndex: line 1 expected 6 columns, got 5"), std::string::npos);
    }
}

TEST(CraiIndex, ThrowsOnNonNumericFieldWithLineAndName)
{
    const TempFileGuard input{
        WriteGzipText("crai_non_numeric", "0\t1\t2\t3\t4\t5\n1\txyz\t2\t3\t4\t5\n")};

    try {
        (void)CraiIndex::FromFile(input.Path());
        FAIL() << "Expected parsing to throw";
    } catch (const std::runtime_error& e) {
        const std::string message = e.what();
        EXPECT_NE(message.find("CraiIndex: line 2 field alignment_start: not numeric: 'xyz'"),
                  std::string::npos);
    }
}

TEST(CraiIndex, ThrowsOnTruncatedGzip)
{
    const std::string line = "0\t1\t2\t3\t4\t5\n";
    const std::span<const std::byte> payload{reinterpret_cast<const std::byte*>(line.data()),
                                             std::size(line)};
    std::vector<std::byte> compressed = CramGzipCompress(payload);
    ASSERT_GT(std::size(compressed), 4U);
    compressed.resize(std::size(compressed) - 4);

    const std::filesystem::path path = TempPath("crai_truncated");
    const TempFileGuard input{path};
    WriteBinary(path, compressed);

    try {
        (void)CraiIndex::FromFile(input.Path());
        FAIL() << "Expected parsing to throw";
    } catch (const std::runtime_error& e) {
        const std::string message = e.what();
        EXPECT_NE(message.find("CraiIndex: failed to decompress gzip data in "), std::string::npos);
        EXPECT_NE(message.find(path.string()), std::string::npos);
        EXPECT_NE(message.find("(truncated or invalid)"), std::string::npos);
    }
}

}  // namespace Samoa
}  // namespace PacBio
