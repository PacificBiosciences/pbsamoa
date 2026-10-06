#include "TestData.hpp"
#include "TestTempDir.hpp"

#include <pbsamoa/core/BamRecord.hpp>
#include <pbsamoa/core/Bgzf.hpp>
#include <pbsamoa/core/SamHeader.hpp>

#include <pbsamoa/io/BamRawReader.hpp>
#include <pbsamoa/io/BamWriter.hpp>

#include <gtest/gtest.h>

#include <algorithm>
#include <filesystem>
#include <format>
#include <fstream>
#include <limits>
#include <optional>
#include <ostream>
#include <span>
#include <string>
#include <vector>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {
namespace {

constexpr std::size_t NO_LIMIT{std::numeric_limits<std::size_t>::max()};

// Writes one unmapped record for each sequence length.
void WriteUnmappedBam(const std::filesystem::path& path, const std::vector<std::size_t>& lengths)
{
    SamHeader header;
    header.SetVersion("1.6");
    BamWriter writer{path, header};
    for (std::size_t i{0}; i < std::size(lengths); ++i) {
        const std::size_t length{lengths[i]};
        BamRecord record;
        record.Name(std::format("r{}", i))
            .Flag(0x4)
            .RefId(-1)
            .Pos(-1)
            .NextRefId(-1)
            .NextPos(-1)
            .Sequence(std::string(length, 'C'))
            .Qualities(std::vector<std::uint8_t>(length, 20));
        writer.Write(record);
    }
}

// The uncompressed BAM stream of a file: the header and the records.
std::vector<std::byte> ReadStream(const std::filesystem::path& path)
{
    BgzfReader reader{path};
    std::vector<std::byte> block(65536);
    std::vector<std::byte> stream;
    for (std::optional<std::size_t> n{reader.ReadBlock(block)}; n && (*n > 0);
         n = reader.ReadBlock(block)) {
        stream.insert(std::end(stream), std::begin(block), std::begin(block) + *n);
    }
    return stream;
}

// Writes the stream as BGZF blocks of blockSize bytes. With tiny blocks, the group
// boundaries of the parallel reader fall at all positions in a record, also in its
// block_size field.
void WriteBgzf(const std::filesystem::path& path, const std::span<const std::byte> stream,
               const std::size_t blockSize)
{
    std::ofstream out{path, std::ios::binary};
    for (std::size_t pos{0}; pos < std::size(stream); pos += blockSize) {
        const std::vector<std::byte> block{CompressBgzfBlock(
            stream.subspan(pos, std::min(blockSize, std::size(stream) - pos)), 1)};
        out.write(reinterpret_cast<const char*>(std::data(block)), std::ssize(block));
    }
    out.write(reinterpret_cast<const char*>(std::data(BGZF_EOF_MARKER)),
              std::ssize(BGZF_EOF_MARKER));
}

struct RecordStream
{
    std::vector<std::size_t> Sizes;
    std::vector<std::byte> Data;

    bool operator==(const RecordStream&) const = default;
};

void PrintTo(const RecordStream& records, std::ostream* os)
{
    *os << std::format("{} records, {} bytes", std::size(records.Sizes), std::size(records.Data));
}

void Append(RecordStream& out, const std::span<const std::byte> record)
{
    out.Sizes.push_back(std::size(record));
    out.Data.insert(std::end(out.Data), std::begin(record), std::end(record));
}

void Append(RecordStream& out, const RawRecordBatch& batch)
{
    for (std::size_t i{0}; i < batch.RecordCount(); ++i) {
        Append(out, batch.RecordData(i));
    }
}

RecordStream ReadAllBatches(BgzfReader& reader, const std::size_t maxBytes)
{
    RecordStream out;
    while (const std::optional<RawRecordBatch> batch{reader.ReadBatch(maxBytes, NO_LIMIT)}) {
        Append(out, *batch);
    }
    return out;
}

RecordStream ReadAllBatches(BamRawReader& reader)
{
    RecordStream out;
    while (const std::optional<RawRecordBatch> batch{reader.ReadBatch()}) {
        Append(out, *batch);
    }
    return out;
}

}  // namespace

TEST(BgzfReaderBam, OpenValidBam)
{
    const std::filesystem::path path{tests::DataDir / "header_only.bam"};
    const BgzfReader pipeline{path, 2};
    EXPECT_TRUE(pipeline.HasEofMarker());
}

TEST(BgzfReaderBam, ThrowOnNonexistentFile)
{
    EXPECT_THROW(BgzfReader("/nonexistent/file.bam", 2), std::runtime_error);
}

TEST(BgzfReaderBam, SyncModeMatchesBgzfReader)
{
    // Compare from file start: BAM-mode reader parses header eagerly.
    const std::filesystem::path path{tests::DataDir / "spec_example.bam"};
    BgzfReader syncReader{path, 0};
    BgzfReader blockReader{path};
    syncReader.Seek(VirtualOffset{0U, 0U});

    std::vector<std::byte> pBuf(65536);
    std::vector<std::byte> rBuf(65536);

    while (true) {
        const auto pResult{syncReader.ReadBlock(pBuf)};
        const auto rResult{blockReader.ReadBlock(rBuf)};

        ASSERT_EQ(pResult.has_value(), rResult.has_value());
        if (!pResult) {
            break;
        }
        ASSERT_EQ(*pResult, *rResult);
        if (*pResult == 0) {
            break;
        }
        EXPECT_TRUE(std::ranges::equal(std::span<const std::byte>{pBuf}.first(*pResult),
                                       std::span<const std::byte>{rBuf}.first(*rResult)));
    }
}

TEST(BgzfReaderBam, PipelineModeMatchesSyncMode)
{
    // Pipeline with N workers must produce identical block-by-block output
    const std::filesystem::path path{tests::DataDir / "spec_example.bam"};
    BgzfReader sync{path, 0};
    BgzfReader parallel{path, 4};

    std::vector<std::byte> sBuf(65536);
    std::vector<std::byte> pBuf(65536);
    std::size_t blockIdx{0};

    while (true) {
        const auto sResult{sync.ReadBlock(sBuf)};
        const auto pResult{parallel.ReadBlock(pBuf)};

        ASSERT_EQ(sResult.has_value(), pResult.has_value()) << "Mismatch at block " << blockIdx;
        if (!sResult) {
            break;
        }
        ASSERT_EQ(*sResult, *pResult) << "Size mismatch at block " << blockIdx;
        if (*sResult == 0) {
            break;
        }
        EXPECT_TRUE(std::ranges::equal(std::span<const std::byte>{sBuf}.first(*sResult),
                                       std::span<const std::byte>{pBuf}.first(*pResult)))
            << "Content mismatch at block " << blockIdx;
        ++blockIdx;
    }
    EXPECT_GT(blockIdx, 0U);
}

TEST(BgzfReaderBam, ManyRecordsBamCorrectness)
{
    // Larger file with many BGZF blocks
    const std::filesystem::path path{tests::DataDir / "many_records.bam"};
    BgzfReader sync{path, 0};
    BgzfReader parallel{path, 8};

    std::vector<std::byte> sBuf(65536);
    std::vector<std::byte> pBuf(65536);
    std::size_t blockCount{0};

    while (true) {
        const auto sResult{sync.ReadBlock(sBuf)};
        const auto pResult{parallel.ReadBlock(pBuf)};

        ASSERT_EQ(sResult.has_value(), pResult.has_value());
        if (!sResult || *sResult == 0) {
            break;
        }
        ASSERT_EQ(*sResult, *pResult);
        EXPECT_TRUE(std::ranges::equal(std::span<const std::byte>{sBuf}.first(*sResult),
                                       std::span<const std::byte>{pBuf}.first(*pResult)));
        ++blockCount;
    }
    EXPECT_GT(blockCount, 1U);
}

TEST(BgzfReaderBam, EofOnSmallFile)
{
    // BAM-mode reader parses header in ctor, so next block is EOF.
    const std::filesystem::path path{tests::DataDir / "header_only.bam"};
    BgzfReader pipeline{path, 2};

    std::vector<std::byte> buf(65536);
    const auto first{pipeline.ReadBlock(buf)};
    ASSERT_TRUE(first);
    EXPECT_EQ(*first, 0U);
}

TEST(BgzfReaderBam, SeekBackToStart)
{
    const std::filesystem::path path{tests::DataDir / "spec_example.bam"};
    BgzfReader pipeline{path, 4};  // starts after header
    BgzfReader blockReader{path};

    std::vector<std::byte> expectedBuf(65536);
    const auto expectedResult{blockReader.ReadBlock(expectedBuf)};
    ASSERT_TRUE(expectedResult);
    ASSERT_GT(*expectedResult, 0U);

    std::vector<std::byte> firstBuf(65536);
    const auto firstResult{pipeline.ReadBlock(firstBuf)};
    ASSERT_TRUE(firstResult);
    ASSERT_GT(*firstResult, 0U);

    // Seek back to start
    pipeline.Seek(VirtualOffset{0U, 0U});

    std::vector<std::byte> secondBuf(65536);
    const auto secondResult{pipeline.ReadBlock(secondBuf)};
    ASSERT_TRUE(secondResult);
    ASSERT_EQ(*expectedResult, *secondResult);
    EXPECT_TRUE(std::ranges::equal(std::span<const std::byte>{expectedBuf}.first(*expectedResult),
                                   std::span<const std::byte>{secondBuf}.first(*secondResult)));
}

TEST(BgzfReaderBam, TellTracksPosition)
{
    const std::filesystem::path path{tests::DataDir / "spec_example.bam"};
    BgzfReader pipeline{path, 2};

    // Before first read, tell should be 0
    EXPECT_EQ(pipeline.Tell().BlockOffset(), 0U);

    std::vector<std::byte> buf(65536);
    const auto result{pipeline.ReadBlock(buf)};
    ASSERT_TRUE(result);
    ASSERT_GT(*result, 0U);

    // After reading first block, tell should be at block 0 offset
    // (consumer reports offset of last consumed block)
    const std::uint64_t afterFirst{pipeline.Tell().BlockOffset()};

    // Read another block
    pipeline.ReadBlock(buf);
    // Tell should have advanced (or be at EOF)
    // For this file we just verify it doesn't regress
    EXPECT_GE(pipeline.Tell().BlockOffset(), afterFirst);
}

TEST(BgzfReaderBam, SingleWorkerCorrectness)
{
    const std::filesystem::path path{tests::DataDir / "spec_example.bam"};
    BgzfReader one{path, 1};
    BgzfReader many{path, 8};

    std::vector<std::byte> oneBuf(65536);
    std::vector<std::byte> manyBuf(65536);

    while (true) {
        const auto oneR{one.ReadBlock(oneBuf)};
        const auto manyR{many.ReadBlock(manyBuf)};

        ASSERT_EQ(oneR.has_value(), manyR.has_value());
        if (!oneR || *oneR == 0) {
            break;
        }
        ASSERT_EQ(*oneR, *manyR);
        EXPECT_TRUE(std::ranges::equal(std::span<const std::byte>{oneBuf}.first(*oneR),
                                       std::span<const std::byte>{manyBuf}.first(*manyR)));
    }
}

// --- ReadRecord tests (v2 pipeline) ---

TEST(BgzfReaderBam, ReadRecordMatchesSyncBamRawReader)
{
    // Pipeline ReadRecord() must produce identical records to sync BamRawReader
    const std::filesystem::path path{tests::DataDir / "spec_example.bam"};

    // Reference: sync BamRawReader reads all records
    BamRawReader syncReader{path};
    std::vector<std::vector<std::byte>> refRecords;
    while (const auto view = syncReader.ReadRecord()) {
        refRecords.emplace_back(std::ranges::begin(view->RawData()),
                                std::ranges::end(view->RawData()));
    }
    ASSERT_FALSE(refRecords.empty());

    // Pipeline: reads records via new API
    BgzfReader pipeline{path, 4};
    std::size_t idx{0};
    while (const auto rec = pipeline.ReadRecord()) {
        ASSERT_LT(idx, refRecords.size()) << "Pipeline produced more records than sync reader";
        EXPECT_TRUE(std::ranges::equal(rec->RawData(), std::span<const std::byte>{refRecords[idx]}))
            << "Record " << idx << " differs";
        ++idx;
    }
    EXPECT_EQ(idx, refRecords.size()) << "Pipeline produced fewer records";
}

TEST(BgzfReaderBam, ReadRecordManyRecordsBam)
{
    const std::filesystem::path path{tests::DataDir / "many_records.bam"};

    BamRawReader syncReader{path};
    std::vector<std::vector<std::byte>> refRecords;
    while (const auto view = syncReader.ReadRecord()) {
        refRecords.emplace_back(std::ranges::begin(view->RawData()),
                                std::ranges::end(view->RawData()));
    }
    ASSERT_FALSE(refRecords.empty());

    BgzfReader pipeline{path, 4};
    std::size_t idx{0};
    while (const auto rec = pipeline.ReadRecord()) {
        ASSERT_LT(idx, refRecords.size());
        EXPECT_TRUE(
            std::ranges::equal(rec->RawData(), std::span<const std::byte>{refRecords[idx]}));
        ++idx;
    }
    EXPECT_EQ(idx, refRecords.size());
}

TEST(BgzfReaderBam, ReadRecordSingleWorker)
{
    const std::filesystem::path path{tests::DataDir / "spec_example.bam"};

    BamRawReader syncReader{path};
    std::vector<std::vector<std::byte>> refRecords;
    while (const auto view = syncReader.ReadRecord()) {
        refRecords.emplace_back(std::ranges::begin(view->RawData()),
                                std::ranges::end(view->RawData()));
    }
    ASSERT_FALSE(refRecords.empty());

    BgzfReader pipeline{path, 1};
    std::size_t idx{0};
    while (const auto rec = pipeline.ReadRecord()) {
        ASSERT_LT(idx, refRecords.size());
        EXPECT_TRUE(
            std::ranges::equal(rec->RawData(), std::span<const std::byte>{refRecords[idx]}));
        ++idx;
    }
    EXPECT_EQ(idx, refRecords.size());
}

TEST(BgzfReaderBam, ReadRecordEofOnHeaderOnly)
{
    const std::filesystem::path path{tests::DataDir / "header_only.bam"};
    BgzfReader pipeline{path, 2};

    // Should immediately return nullopt — no records in header-only file
    const auto rec{pipeline.ReadRecord()};
    EXPECT_FALSE(rec);
}

TEST(BgzfReaderBam, ParallelReadSurfacesTruncatedRecordBlock)
{
    // A truncated BGZF block in the record region must surface as a read error in
    // parallel mode, not finalize as a clean EOF (which would silently drop the records
    // already accumulated in the in-flight batch). Copy a valid multi-block BAM and cut
    // into the last data block by dropping the 28-byte EOF marker plus two body bytes.
    const std::filesystem::path src{tests::DataDir / "many_records.bam"};
    tests::TempDirGuard tempDir{"bgzf_trunc"};
    const std::filesystem::path corrupt{tempDir.File("truncated.bam")};
    std::filesystem::copy_file(src, corrupt);
    const auto size{std::filesystem::file_size(corrupt)};
    ASSERT_GT(size, 30U);
    std::filesystem::resize_file(corrupt, size - 30);

    BgzfReader pipeline{corrupt, 4};
    EXPECT_THROW(
        {
            while (pipeline.ReadRecord()) {
            }
        },
        std::runtime_error);
}

TEST(BgzfReaderBam, ParallelReadBatchMatchesSyncAtAllGroupBoundaries)
{
    const tests::TempDirGuard tempDir{"bgzf_tiny_blocks"};
    const std::filesystem::path source{tempDir.File("source.bam")};
    // The short records put a group boundary in a block_size field 24 times, most of
    // them with blocks of 1 and 3 bytes.
    std::vector<std::size_t> lengths{1, 2, 3, 10, 50, 333, 1000, 2000, 7, 1, 1500, 4, 900};
    for (std::size_t length{1}; length <= 64; ++length) {
        lengths.push_back(length);
    }
    WriteUnmappedBam(source, lengths);
    const std::vector<std::byte> stream{ReadStream(source)};

    for (const std::size_t blockSize : {1, 3, 7, 64, 1000}) {
        SCOPED_TRACE(blockSize);
        const std::filesystem::path path{tempDir.File(std::format("blocks{}.bam", blockSize))};
        WriteBgzf(path, stream, blockSize);

        BgzfReader sync{path, 0};
        const RecordStream expected{ReadAllBatches(sync, NO_LIMIT)};
        ASSERT_EQ(std::size(expected.Sizes), std::size(lengths));

        // The batch takes over the group buffer.
        BgzfReader parallel{path, 4};
        EXPECT_EQ(ReadAllBatches(parallel, NO_LIMIT), expected);

        // One record for each batch: the records are copied out of the group.
        BgzfReader parallelOneRecord{path, 4};
        EXPECT_EQ(ReadAllBatches(parallelOneRecord, 1), expected);
    }
}

TEST(BgzfReaderBam, ParallelReadBatchKeepsARecordLargerThanAGroup)
{
    // A group is 32 blocks, about 2 MiB. The 2 Mbp record (about 3 MB) spans two groups
    // and is larger than the spare capacity of the group where it ends.
    const tests::TempDirGuard tempDir{"bgzf_large_record"};
    const std::filesystem::path path{tempDir.File("large.bam")};
    WriteUnmappedBam(path, {10, 2'000'000, 10, 100'000, 10});

    BgzfReader sync{path, 0};
    const RecordStream expected{ReadAllBatches(sync, NO_LIMIT)};
    ASSERT_EQ(std::size(expected.Sizes), 5U);

    BgzfReader parallel{path, 4};
    EXPECT_EQ(ReadAllBatches(parallel, NO_LIMIT), expected);
}

TEST(BgzfReaderBam, ParallelReadRecordAndReadBatchMix)
{
    const tests::TempDirGuard tempDir{"bgzf_mix"};
    const std::filesystem::path path{tempDir.File("mix.bam")};
    WriteBgzf(path, ReadStream(tests::DataDir / "many_records.bam"), 512);

    BgzfReader sync{path, 0};
    const RecordStream expected{ReadAllBatches(sync, NO_LIMIT)};
    ASSERT_EQ(std::size(expected.Sizes), 500U);

    // A record, then a batch of at most 2 records, then the rest of the group.
    BgzfReader parallel{path, 4};
    RecordStream actual;
    for (std::size_t step{0};; ++step) {
        if ((step % 3) == 0) {
            const std::optional<RawRecord> record{parallel.ReadRecord()};
            if (!record) {
                break;
            }
            Append(actual, record->RawData());
            continue;
        }
        const std::size_t maxRecords{((step % 3) == 1) ? std::size_t{2} : NO_LIMIT};
        const std::optional<RawRecordBatch> batch{parallel.ReadBatch(NO_LIMIT, maxRecords)};
        if (!batch) {
            break;
        }
        EXPECT_LE(batch->RecordCount(), maxRecords);
        Append(actual, *batch);
    }
    EXPECT_EQ(actual, expected);
}

TEST(BgzfReaderBam, ParallelRecordLimitCopiesPartOfAGroup)
{
    // A group of 32 blocks of 512 bytes holds about 46 records, so the record limit
    // takes the first 7 records of the first group.
    const tests::TempDirGuard tempDir{"bgzf_record_limit"};
    const std::filesystem::path path{tempDir.File("limit.bam")};
    WriteBgzf(path, ReadStream(tests::DataDir / "many_records.bam"), 512);

    BamRawReader sync{path, BamRawReaderConfig{.RecordLimit = 7}};
    const RecordStream expected{ReadAllBatches(sync)};
    ASSERT_EQ(std::size(expected.Sizes), 7U);

    BamRawReader parallel{path, BamRawReaderConfig{.BgzfWorkers = 4, .RecordLimit = 7}};
    EXPECT_EQ(ReadAllBatches(parallel), expected);
}

TEST(BgzfReaderBam, ParseHeaderMatchesBamRawReader)
{
    const std::filesystem::path path{tests::DataDir / "spec_example.bam"};

    const BamRawReader syncReader{path};
    const BgzfReader pipeline{path, 4};

    EXPECT_EQ(pipeline.Header().NumReferences(), syncReader.Header().NumReferences());
    EXPECT_EQ(pipeline.Header().Version(), syncReader.Header().Version());
}

}  // namespace Samoa
}  // namespace PacBio
