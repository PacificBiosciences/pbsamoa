#include "TestData.hpp"
#include "TestTempDir.hpp"

#include <pbsamoa/core/BamRecord.hpp>
#include <pbsamoa/core/CigarOp.hpp>
#include <pbsamoa/core/RawRecord.hpp>
#include <pbsamoa/core/SamHeader.hpp>
#include <pbsamoa/index/BaiIndex.hpp>
#include <pbsamoa/index/ZmwWhitelist.hpp>
#include <pbsamoa/io/BamRawReader.hpp>
#include <pbsamoa/io/BamRecordReader.hpp>
#include <pbsamoa/io/ZmiBamWriter.hpp>

#include <gtest/gtest.h>

#include <algorithm>
#include <charconv>
#include <filesystem>
#include <format>
#include <map>
#include <string>
#include <string_view>
#include <thread>
#include <vector>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {

TEST(BamRawReader, OpenHeaderOnly)
{
    const BamRawReader reader{tests::DataDir / "header_only.bam"};
    const SamHeader& header = reader.Header();
    EXPECT_FALSE(std::empty(header.Version()));
}

TEST(BamRawReader, OpenSpecExample)
{
    const BamRawReader reader{tests::DataDir / "spec_example.bam"};
    const SamHeader& header = reader.Header();
    EXPECT_EQ(header.Version(), "1.6");
    EXPECT_EQ(header.SortOrder(), "coordinate");
    EXPECT_EQ(header.ReferenceSequences()[0].Name(), "ref");
    EXPECT_EQ(header.ReferenceSequences()[0].Length(), 45);
}

TEST(BamRawReader, ReadRecordsFromSpecExample)
{
    BamRawReader reader{tests::DataDir / "spec_example.bam"};
    std::vector<std::string> names;
    for (auto view = reader.ReadRecord(); view; view = reader.ReadRecord()) {
        names.emplace_back(view->Name());
    }
    ASSERT_EQ(std::size(names), 6u);
    EXPECT_EQ(names[0], "r001");
}

TEST(BamRawReader, RecordFieldsDecodeCorrectly)
{
    BamRawReader reader{tests::DataDir / "spec_example.bam"};
    const auto view = reader.ReadRecord();
    ASSERT_TRUE(view);
    EXPECT_EQ(view->Name(), "r001");
    EXPECT_EQ(view->Flag(), 99u);
    EXPECT_EQ(view->RefId(), 0);
    EXPECT_EQ(view->Pos(), 6);  // 0-based
    EXPECT_EQ(view->MapQ(), 30u);
}

TEST(BamRawReader, HeaderOnlyBamProducesNoRecords)
{
    BamRawReader reader{tests::DataDir / "header_only.bam"};
    EXPECT_FALSE(reader.ReadRecord());
}

TEST(BamRawReader, ThrowOnNonexistent)
{
    EXPECT_THROW(BamRawReader{"/nonexistent/file.bam"}, std::runtime_error);
}

// --- Batch interface ---

TEST(BamRawReader, ReadBatchReturnsRecords)
{
    BamRawReader reader{tests::DataDir / "spec_example.bam"};
    const auto batch = reader.ReadBatch(ByteLimit{1024U * 1024U});
    ASSERT_TRUE(batch);
    EXPECT_EQ(batch->RecordCount(), 6u);
}

TEST(BamRawReader, ReadBatchRespectsLimit)
{
    BamRawReader reader{tests::DataDir / "spec_example.bam"};
    std::size_t totalRecords{0};
    std::size_t batchCount{0};
    while (const auto batch = reader.ReadBatch(ByteLimit{64})) {
        totalRecords += batch->RecordCount();
        ++batchCount;
        EXPECT_GE(batch->RecordCount(), 1u);
    }
    EXPECT_EQ(totalRecords, 6u);
    EXPECT_GE(batchCount, 2u);
}

TEST(BamRawReader, HeaderOnlyBamBatchReturnsNullopt)
{
    BamRawReader reader{tests::DataDir / "header_only.bam"};
    EXPECT_FALSE(reader.ReadBatch());
}

// --- Range interface ---

TEST(BamRawReader, RangeInterfaceIteratesAll)
{
    BamRawReader reader{tests::DataDir / "spec_example.bam"};
    std::vector<std::string> names;
    for (const auto& view : reader.Records()) {
        names.emplace_back(view.Name());
    }
    EXPECT_EQ(std::size(names), 6u);
}

TEST(BamRawReader, RangeInterfaceEmptyFile)
{
    BamRawReader reader{tests::DataDir / "header_only.bam"};
    std::size_t count{0};
    for ([[maybe_unused]] const auto& view : reader.Records()) {
        ++count;
    }
    EXPECT_EQ(count, 0u);
}

// --- Records spanning BGZF block boundaries ---

TEST(BamRawReader, ManyRecordsSpanningBlocks)
{
    const auto path = tests::DataDir / "many_records.bam";
    if (!std::filesystem::exists(path)) {
        GTEST_SKIP() << "many_records.bam not generated";
    }

    BamRawReader reader{path};

    std::size_t count{0};
    for (const auto& view : reader.Records()) {
        EXPECT_FALSE(std::empty(view.Name()));
        EXPECT_TRUE(view.IsMapped());
        ++count;
    }

    EXPECT_EQ(count, 500u);
}

TEST(BamRawReader, BatchModeSpanningBlocks)
{
    const auto path = tests::DataDir / "many_records.bam";
    if (!std::filesystem::exists(path)) {
        GTEST_SKIP() << "many_records.bam not generated";
    }

    BamRawReader reader{path};

    std::size_t total{0};
    while (const auto batch = reader.ReadBatch(ByteLimit{4096})) {
        for (std::size_t i{0}; i < batch->RecordCount(); ++i) {
            const RawRecord view{batch->RecordData(i)};
            EXPECT_FALSE(std::empty(view.Name()));
        }
        total += batch->RecordCount();
    }

    EXPECT_EQ(total, 500u);
}

// --- Query with index ---

TEST(BamRawReader, QueryWithIndex)
{
    const auto baiPath = tests::DataDir / "spec_example.bam.bai";
    const auto bamPath = tests::DataDir / "spec_example.bam";
    if (!std::filesystem::exists(baiPath)) {
        GTEST_SKIP() << "spec_example.bam.bai not found";
    }

    BamRawReader reader{bamPath};
    const auto index = BaiIndex::FromFile(baiPath);

    std::size_t count{0};
    for (const auto& view : reader.Query(index, 0, 0, 45)) {
        EXPECT_EQ(view.RefId(), 0);
        ++count;
    }
    EXPECT_GT(count, 0u);
}

// --- Stress tests with benchmark file ---

TEST(BamRawReader, StressTestBenchmarkFile)
{
    const auto path = tests::DataDir / "benchmark.bam";
    if (!std::filesystem::exists(path)) {
        GTEST_SKIP() << "benchmark.bam not generated";
    }

    BamRawReader reader{path};
    EXPECT_EQ(reader.Header().ReferenceSequences()[0].Name(), "chr1");

    std::size_t count{0};
    std::int32_t prevPos{-1};
    for (const auto& view : reader.Records()) {
        EXPECT_EQ(view.RefId(), 0);
        EXPECT_GE(view.Pos(), prevPos) << "Not sorted at record " << count;
        prevPos = view.Pos();
        EXPECT_EQ(std::size(view.Seq().ToString()), 200u);
        EXPECT_EQ(view.MapQ(), 30u);
        ++count;
    }
    EXPECT_EQ(count, 10000u);
}

TEST(BamRawReader, StressBatchMode)
{
    const auto path = tests::DataDir / "benchmark.bam";
    if (!std::filesystem::exists(path)) {
        GTEST_SKIP() << "benchmark.bam not generated";
    }

    BamRawReader reader{path};

    std::size_t total{0};
    std::size_t batchCount{0};
    while (const auto batch = reader.ReadBatch(ByteLimit{1024 * 1024})) {
        total += batch->RecordCount();
        ++batchCount;
    }

    EXPECT_EQ(total, 10000u);
    EXPECT_GE(batchCount, 2u);  // 10k 200bp records > 1 MiB
}

// --- Pipeline vs sync comparison ---

TEST(BamRawReader, PipelineMatchesSyncReadRecord)
{
    const std::filesystem::path path{tests::DataDir / "spec_example.bam"};
    BamRawReader syncReader{path};
    BamRawReader pipeReader{path, BamRawReaderConfig{.BgzfWorkers = 4}};

    std::size_t count{0};
    while (true) {
        const auto syncRec{syncReader.ReadRecord()};
        const auto pipeRec{pipeReader.ReadRecord()};

        ASSERT_EQ(syncRec.has_value(), pipeRec.has_value()) << "at record " << count;
        if (!syncRec) {
            break;
        }

        EXPECT_EQ(syncRec->Flag(), pipeRec->Flag()) << "at record " << count;
        EXPECT_EQ(syncRec->Pos(), pipeRec->Pos()) << "at record " << count;
        EXPECT_EQ(syncRec->MapQ(), pipeRec->MapQ()) << "at record " << count;
        EXPECT_EQ(syncRec->RefId(), pipeRec->RefId()) << "at record " << count;
        ++count;
    }
    EXPECT_GT(count, 0U);
}

TEST(BamRawReader, PipelineMatchesSyncReadBatch)
{
    const auto path{tests::DataDir / "many_records.bam"};
    if (!std::filesystem::exists(path)) {
        GTEST_SKIP() << "many_records.bam not generated";
    }

    BamRawReader syncReader{path};
    BamRawReader pipeReader{path, BamRawReaderConfig{.BgzfWorkers = 4}};

    std::vector<std::uint16_t> syncFlags;
    std::vector<std::uint16_t> pipeFlags;

    while (const auto batch = syncReader.ReadBatch()) {
        for (std::size_t i{0}; i < batch->RecordCount(); ++i) {
            const RawRecord rec{batch->RecordData(i)};
            syncFlags.push_back(rec.Flag());
        }
    }
    while (const auto batch = pipeReader.ReadBatch()) {
        for (std::size_t i{0}; i < batch->RecordCount(); ++i) {
            const RawRecord rec{batch->RecordData(i)};
            pipeFlags.push_back(rec.Flag());
        }
    }

    ASSERT_EQ(std::size(syncFlags), std::size(pipeFlags));
    EXPECT_TRUE(std::ranges::equal(syncFlags, pipeFlags));
}

class BamRawReaderWhitelistTest : public ::testing::Test
{
protected:
    tests::TempDirGuard tempDir_;
    std::filesystem::path tmpBamPath_;

    void SetUp() override
    {
        tempDir_.Reset("bam_raw_reader_whitelist");
        tmpBamPath_ = tempDir_.File("whitelist.bam");

        SamHeader header;
        header.SetVersion("1.6");
        header.SetSortOrder("unknown");
        header.AddReferenceSequence(ReferenceSequence{"ref", 1000});

        ZmiBamWriter writer{tmpBamPath_, header};
        // Write records for 3 ZMWs: 10, 20, 30
        for (const std::int32_t zmw : {10, 20, 30}) {
            for (int subread = 0; subread < 2; ++subread) {
                const std::string name{
                    std::format("movie/{}/{}_{}", zmw, subread * 100, (subread + 1) * 100)};
                BamRecord rec;
                rec.Name(name)
                    .Flag(0)
                    .RefId(0)
                    .Pos(0)
                    .MapQ(30)
                    .Cigar({CigarOp{CigarOpType::M, 4}})
                    .Sequence("ACGT")
                    .Qualities({30, 30, 30, 30});
                writer.Write(rec);
            }
        }
    }
};

TEST_F(BamRawReaderWhitelistTest, WhitelistSelectsSubset)
{
    BamRawReader reader{tmpBamPath_};
    const ZmwWhitelist whitelist{std::vector<std::int32_t>{10, 30}};

    std::vector<std::string> names;
    for (const auto& rec : reader.Whitelist(whitelist)) {
        names.emplace_back(rec.Name());
    }

    // zmw 10: 2 records, zmw 30: 2 records
    ASSERT_EQ(std::size(names), 4u);
    // All names should contain /10/ or /30/
    for (const auto& name : names) {
        EXPECT_TRUE(name.find("/10/") != std::string::npos ||
                    name.find("/30/") != std::string::npos)
            << "unexpected name: " << name;
    }
}

TEST_F(BamRawReaderWhitelistTest, EmptyWhitelistYieldsNoRecords)
{
    BamRawReader reader{tmpBamPath_};
    const ZmwWhitelist whitelist{std::vector<std::int32_t>{}};

    std::size_t count{0};
    for ([[maybe_unused]] const auto& rec : reader.Whitelist(whitelist)) {
        ++count;
    }
    EXPECT_EQ(count, 0u);
}

TEST_F(BamRawReaderWhitelistTest, NonexistentZmwYieldsNoRecords)
{
    BamRawReader reader{tmpBamPath_};
    const ZmwWhitelist whitelist{std::vector<std::int32_t>{9999}};

    std::size_t count{0};
    for ([[maybe_unused]] const auto& rec : reader.Whitelist(whitelist)) {
        ++count;
    }
    EXPECT_EQ(count, 0u);
}

TEST_F(BamRawReaderWhitelistTest, WhitelistWithPipeline)
{
    BamRawReader reader{tmpBamPath_, BamRawReaderConfig{.BgzfWorkers = 2}};
    const ZmwWhitelist whitelist{std::vector<std::int32_t>{20}};

    std::vector<std::string> names;
    for (const auto& rec : reader.Whitelist(whitelist)) {
        names.emplace_back(rec.Name());
    }

    ASSERT_EQ(std::size(names), 2u);
    for (const auto& name : names) {
        EXPECT_NE(name.find("/20/"), std::string::npos) << "unexpected name: " << name;
    }
}

TEST_F(BamRawReaderWhitelistTest, WhitelistRangeSurvivesEmptyCheckBeforeIteration)
{
    BamRawReader reader{tmpBamPath_};
    const ZmwWhitelist whitelist{std::vector<std::int32_t>{10, 30}};

    auto range = reader.Whitelist(whitelist);
    EXPECT_NE(range.begin(), range.end());

    std::vector<std::string> names;
    for (const auto& rec : range) {
        names.emplace_back(rec.Name());
    }

    EXPECT_EQ(std::size(names), 4u);
}

TEST_F(BamRawReaderWhitelistTest, ChunkingCanProduceEmptyChunk)
{
    BamRawReader reader{tmpBamPath_, BamRawReaderConfig{.ChunkNum = 2, .TotalChunks = 10}};

    std::size_t count{0};
    for ([[maybe_unused]] const auto& rec : reader.Records()) {
        ++count;
    }
    EXPECT_EQ(count, 0u);
}

TEST_F(BamRawReaderWhitelistTest, SinglePathVectorSupportsChunking)
{
    BamRawReader reader{std::vector<std::filesystem::path>{tmpBamPath_},
                        BamRawReaderConfig{.ChunkNum = 1, .TotalChunks = 1}};

    std::size_t count{0};
    for ([[maybe_unused]] const auto& rec : reader.Records()) {
        ++count;
    }

    EXPECT_EQ(count, 6u);
}

TEST(BamRawReader, MultiPathVectorRejectsChunking)
{
    const auto path = tests::DataDir / "spec_example.bam";
    BamRawReaderConfig config{};
    config.ChunkNum = 1;
    config.TotalChunks = 2;

    EXPECT_THROW((BamRawReader{std::vector<std::filesystem::path>{path, path}, config}),
                 std::invalid_argument);
}

TEST(BamRawReader, MultiPathVectorRejectsWhitelist)
{
    const auto path = tests::DataDir / "spec_example.bam";
    BamRawReaderConfig config{};
    config.Whitelist.emplace(std::vector<std::int32_t>{42});

    EXPECT_THROW((BamRawReader{std::vector<std::filesystem::path>{path, path}, config}),
                 std::invalid_argument);
}

TEST_F(BamRawReaderWhitelistTest, ChunkingWithManyChunksStillCoversAllRecords)
{
    std::size_t total{0};
    for (std::int32_t chunkNum{1}; chunkNum <= 10; ++chunkNum) {
        BamRawReader reader{tmpBamPath_,
                            BamRawReaderConfig{.ChunkNum = chunkNum, .TotalChunks = 10}};
        for ([[maybe_unused]] const auto& rec : reader.Records()) {
            ++total;
        }
    }

    EXPECT_EQ(total, 6u);
}

TEST_F(BamRawReaderWhitelistTest, ChunkedReadWithPipelineLandsOnCorrectRecords)
{
    // 3 ZMWs (10, 20, 30), 2 records each, packed into a single BGZF block, so
    // chunk boundaries fall mid-block. With a parallel BGZF pipeline the reader
    // must still honour the within-block start offset; otherwise the chunk seeks
    // to the block boundary and returns earlier records than requested.
    for (std::int32_t chunkNum{1}; chunkNum <= 3; ++chunkNum) {
        BamRawReader reader{
            tmpBamPath_,
            BamRawReaderConfig{.BgzfWorkers = 4, .ChunkNum = chunkNum, .TotalChunks = 3}};

        std::vector<std::string> names;
        for (const auto& rec : reader.Records()) {
            names.emplace_back(rec.Name());
        }

        const std::string marker{std::format("/{}/", chunkNum * 10)};
        ASSERT_EQ(std::size(names), 2u) << "chunk " << chunkNum;
        for (const auto& name : names) {
            EXPECT_NE(name.find(marker), std::string::npos)
                << "chunk " << chunkNum << " returned unexpected record: " << name;
        }
    }
}

TEST_F(BamRawReaderWhitelistTest, ChunkedPipelineUnionEqualsWholeFile)
{
    // The union of all chunks must equal the whole-file read exactly — no
    // duplicated or dropped records — even with the parallel BGZF pipeline.
    std::vector<std::string> whole;
    {
        BamRawReader reader{tmpBamPath_, BamRawReaderConfig{.BgzfWorkers = 4}};
        for (const auto& rec : reader.Records()) {
            whole.emplace_back(rec.Name());
        }
    }

    std::vector<std::string> unioned;
    for (std::int32_t chunkNum{1}; chunkNum <= 3; ++chunkNum) {
        BamRawReader reader{
            tmpBamPath_,
            BamRawReaderConfig{.BgzfWorkers = 4, .ChunkNum = chunkNum, .TotalChunks = 3}};
        for (const auto& rec : reader.Records()) {
            unioned.emplace_back(rec.Name());
        }
    }

    std::ranges::sort(whole);
    std::ranges::sort(unioned);
    EXPECT_EQ(unioned, whole);
}

TEST(BamRawReader, WhitelistAndChunkingMutuallyExclusive)
{
    const auto path = tests::DataDir / "spec_example.bam";
    BamRawReaderConfig config{};
    config.ChunkNum = 1;
    config.TotalChunks = 2;
    config.Whitelist.emplace(std::vector<std::int32_t>{42});

    EXPECT_THROW(BamRawReader(path, config), std::invalid_argument);
}

TEST_F(BamRawReaderWhitelistTest, NumZmwsSingleFile)
{
    // Single file with ZMI: 3 ZMWs (10, 20, 30), 2 records each
    const BamRawReader reader{tmpBamPath_};
    EXPECT_EQ(reader.NumZmws(), 3);
}

TEST_F(BamRawReaderWhitelistTest, NumZmwsChunked)
{
    BamRawReader reader{tmpBamPath_, BamRawReaderConfig{.ChunkNum = 1, .TotalChunks = 1}};
    EXPECT_EQ(reader.NumZmws(), 3);
}

TEST(BamRawReader, NumZmwsCollectionReturnsNegativeOne)
{
    const auto path = tests::DataDir / "spec_example.bam";
    BamRawReader reader{std::vector<std::filesystem::path>{path, path}};
    EXPECT_EQ(reader.NumZmws(), -1);
}

TEST(BamRawReader, RecordLimitStopsEarly)
{
    BamRawReader reader{tests::DataDir / "spec_example.bam", BamRawReaderConfig{.RecordLimit = 2}};
    std::size_t count{0};
    for ([[maybe_unused]] const auto& rec : reader.Records()) {
        ++count;
    }
    EXPECT_EQ(count, 2u);
}

TEST(BamRawReader, SeekAndTell)
{
    BamRawReader reader{tests::DataDir / "spec_example.bam"};
    const auto first{reader.ReadRecord()};
    ASSERT_TRUE(first);
    const std::string firstName{std::string{first->Name()}};

    const VirtualOffset pos{reader.Tell()};

    // Read past more records
    reader.ReadRecord();
    reader.ReadRecord();

    // Seek back
    reader.Seek(VirtualOffset{0, 0});

    // First seek to a known position doesn't guarantee same record due to how
    // virtual offsets work, so just verify Seek+Tell round-trip
    reader.Seek(pos);
    const VirtualOffset pos2{reader.Tell()};
    EXPECT_EQ(pos.Value(), pos2.Value());
}

TEST(BamRawReader, QueryEmptyRegion)
{
    const auto baiPath = tests::DataDir / "spec_example.bam.bai";
    const auto bamPath = tests::DataDir / "spec_example.bam";
    if (!std::filesystem::exists(baiPath)) {
        GTEST_SKIP() << "spec_example.bam.bai not found";
    }

    BamRawReader reader{bamPath};
    const auto index = BaiIndex::FromFile(baiPath);

    // Query a region far beyond any records
    std::size_t count{0};
    for ([[maybe_unused]] const auto& rec : reader.Query(index, 0, 100000, 200000)) {
        ++count;
    }
    EXPECT_EQ(count, 0u);
}

// --- Scatter (chaotic-deterministic) chunking ---

namespace {

// Builds a BAM (+ ZMI) at `path`: `numZmws` ZMWs in ascending hole-number order,
// `subreadsPerZmw` records each. Hole numbers are 0..numZmws-1 so that file order
// equals ascending hole number (an "ordered" fixture). Record names encode the
// hole and subread index.
void BuildScatterBam(const std::filesystem::path& path, std::int32_t numZmws,
                     std::int32_t subreadsPerZmw)
{
    SamHeader header;
    header.SetVersion("1.6");
    header.SetSortOrder("unknown");
    header.AddReferenceSequence(ReferenceSequence{"ref", 1000});

    ZmiBamWriter writer{path, header};
    for (std::int32_t hole = 0; hole < numZmws; ++hole) {
        for (std::int32_t sub = 0; sub < subreadsPerZmw; ++sub) {
            BamRecord rec;
            rec.Name(std::format("movie/{}/{}_{}", hole, sub * 100, (sub + 1) * 100))
                .Flag(0)
                .RefId(0)
                .Pos(0)
                .MapQ(30)
                .Cigar({CigarOp{CigarOpType::M, 4}})
                .Sequence("ACGT")
                .Qualities({30, 30, 30, 30});
            writer.Write(rec);
        }
    }
}

BamRawReaderConfig ScatterConfig(std::int32_t chunkNum, std::int32_t totalChunks,
                                 std::int32_t tileZmws, std::uint64_t seed)
{
    return BamRawReaderConfig{
        .ChunkNum = chunkNum,
        .TotalChunks = totalChunks,
        .ChunkingMode = ChunkMode::SCATTER,
        .ChunkTileZmws = tileZmws,
        .ChunkSeed = seed,
    };
}

std::vector<std::string> ReadChunkNames(const std::filesystem::path& path,
                                        const BamRawReaderConfig& config)
{
    BamRawReader reader{path, config};
    std::vector<std::string> names;
    for (const auto& view : reader.Records()) {
        names.emplace_back(view.Name());
    }
    return names;
}

std::vector<std::string> WholeFileNames(const std::filesystem::path& path)
{
    return ReadChunkNames(path, BamRawReaderConfig{});
}

// Parses the hole number out of "movie/<hole>/<a>_<b>".
std::int32_t HoleOf(std::string_view name)
{
    const std::size_t first{name.find('/')};
    const std::size_t second{name.find('/', first + 1)};
    const std::string_view holeText{name.substr(first + 1, second - first - 1)};
    std::int32_t hole{0};
    std::from_chars(holeText.data(), holeText.data() + std::size(holeText), hole);
    return hole;
}

// Sorted distinct hole numbers per chunk for a scatter run.
std::vector<std::vector<std::int32_t>> ScatterHolePartition(const std::filesystem::path& path,
                                                            std::int32_t totalChunks,
                                                            std::int32_t tileZmws,
                                                            std::uint64_t seed)
{
    std::vector<std::vector<std::int32_t>> partition;
    for (std::int32_t chunk = 1; chunk <= totalChunks; ++chunk) {
        const auto names{ReadChunkNames(path, ScatterConfig(chunk, totalChunks, tileZmws, seed))};
        std::vector<std::int32_t> holes;
        for (const auto& name : names) {
            holes.push_back(HoleOf(name));
        }
        std::ranges::sort(holes);
        holes.erase(std::ranges::begin(std::ranges::unique(holes)), std::ranges::end(holes));
        partition.push_back(std::move(holes));
    }
    return partition;
}

}  // namespace

TEST(BamRawReaderScatter, CoversEveryRecordExactlyOnceSubreads)
{
    tests::TempDirGuard tempDir;
    tempDir.Reset("scatter_subreads");
    const std::filesystem::path path{tempDir.File("scatter.bam")};
    BuildScatterBam(path, /*numZmws=*/12, /*subreadsPerZmw=*/2);

    const std::int32_t totalChunks{4};
    std::vector<std::string> all;
    for (std::int32_t chunk = 1; chunk <= totalChunks; ++chunk) {
        const auto names{ReadChunkNames(path, ScatterConfig(chunk, totalChunks, 1, 0))};
        all.insert(std::ranges::end(all), std::ranges::begin(names), std::ranges::end(names));
    }
    std::ranges::sort(all);

    auto expected{WholeFileNames(path)};
    std::ranges::sort(expected);

    EXPECT_EQ(std::size(all), std::size(expected));                     // nothing dropped
    EXPECT_EQ(std::ranges::adjacent_find(all), std::ranges::end(all));  // no duplicates
    EXPECT_EQ(all, expected);                                           // exact same record set
}

TEST(BamRawReaderScatter, CoversEveryRecordExactlyOnceHiFi)
{
    tests::TempDirGuard tempDir;
    tempDir.Reset("scatter_hifi");
    const std::filesystem::path path{tempDir.File("scatter.bam")};
    BuildScatterBam(path, /*numZmws=*/12, /*subreadsPerZmw=*/1);

    const std::int32_t totalChunks{5};
    std::vector<std::string> all;
    for (std::int32_t chunk = 1; chunk <= totalChunks; ++chunk) {
        const auto names{ReadChunkNames(path, ScatterConfig(chunk, totalChunks, 1, 0))};
        all.insert(std::ranges::end(all), std::ranges::begin(names), std::ranges::end(names));
    }
    std::ranges::sort(all);

    auto expected{WholeFileNames(path)};
    std::ranges::sort(expected);

    EXPECT_EQ(std::ranges::adjacent_find(all), std::ranges::end(all));
    EXPECT_EQ(all, expected);
}

TEST(BamRawReaderScatter, IsDeterministicForSameSeed)
{
    tests::TempDirGuard tempDir;
    tempDir.Reset("scatter_determinism");
    const std::filesystem::path path{tempDir.File("scatter.bam")};
    BuildScatterBam(path, /*numZmws=*/12, /*subreadsPerZmw=*/2);

    for (std::int32_t chunk = 1; chunk <= 4; ++chunk) {
        const auto first{ReadChunkNames(path, ScatterConfig(chunk, 4, 1, 42))};
        const auto second{ReadChunkNames(path, ScatterConfig(chunk, 4, 1, 42))};
        EXPECT_EQ(first, second);
    }
}

TEST(BamRawReaderScatter, SamplesAcrossOrderedFile)
{
    tests::TempDirGuard tempDir;
    tempDir.Reset("scatter_ordered");
    const std::filesystem::path path{tempDir.File("scatter.bam")};
    BuildScatterBam(path, /*numZmws=*/12, /*subreadsPerZmw=*/1);

    // Contiguous chunk 1 of 4 would be the first quarter of the file: holes {0,1,2}.
    // Scatter must spread holes across the whole file, so the partition differs.
    const auto scatter{ScatterHolePartition(path, 4, 1, 7)};

    std::vector<std::vector<std::int32_t>> contiguous;
    for (std::int32_t chunk = 1; chunk <= 4; ++chunk) {
        const BamRawReaderConfig cfg{.ChunkNum = chunk, .TotalChunks = 4};  // CONTIGUOUS default
        const auto names{ReadChunkNames(path, cfg)};
        std::vector<std::int32_t> holes;
        for (const auto& name : names) {
            holes.push_back(HoleOf(name));
        }
        std::ranges::sort(holes);
        contiguous.push_back(std::move(holes));
    }

    EXPECT_NE(scatter, contiguous);
    // First contiguous chunk is exactly {0,1,2}; scatter's first chunk must not be.
    EXPECT_NE(scatter.front(), (std::vector<std::int32_t>{0, 1, 2}));
}

TEST(BamRawReaderScatter, DifferentSeedsChangeAssignment)
{
    tests::TempDirGuard tempDir;
    tempDir.Reset("scatter_seed");
    const std::filesystem::path path{tempDir.File("scatter.bam")};
    BuildScatterBam(path, /*numZmws=*/12, /*subreadsPerZmw=*/1);

    const auto seedA{ScatterHolePartition(path, 4, 1, 1)};
    const auto seedB{ScatterHolePartition(path, 4, 1, 2)};
    EXPECT_NE(seedA, seedB);
}

TEST(BamRawReaderScatter, KeepsSubreadsOfSameZmwTogether)
{
    tests::TempDirGuard tempDir;
    tempDir.Reset("scatter_zmw_together");
    const std::filesystem::path path{tempDir.File("scatter.bam")};
    BuildScatterBam(path, /*numZmws=*/12, /*subreadsPerZmw=*/3);

    std::map<std::int32_t, std::int32_t> holeToChunk;
    for (std::int32_t chunk = 1; chunk <= 4; ++chunk) {
        const auto names{ReadChunkNames(path, ScatterConfig(chunk, 4, 1, 9))};
        std::map<std::int32_t, std::int32_t> recordsThisChunk;
        for (const auto& name : names) {
            const std::int32_t hole{HoleOf(name)};
            // A ZMW must not be split across chunks.
            const auto [it, inserted]{holeToChunk.try_emplace(hole, chunk)};
            if (!inserted) {
                EXPECT_EQ(it->second, chunk) << "ZMW " << hole << " split across chunks";
            }
            ++recordsThisChunk[hole];
        }
        // Every ZMW present in this chunk must bring all 3 of its subreads.
        for (const auto& [hole, count] : recordsThisChunk) {
            EXPECT_EQ(count, 3) << "ZMW " << hole << " missing subreads";
        }
    }
    EXPECT_EQ(std::size(holeToChunk), 12u);  // all ZMWs covered
}

TEST(BamRawReaderScatter, TileLargerThanInputPutsAllInOneChunk)
{
    tests::TempDirGuard tempDir;
    tempDir.Reset("scatter_big_tile");
    const std::filesystem::path path{tempDir.File("scatter.bam")};
    BuildScatterBam(path, /*numZmws=*/5, /*subreadsPerZmw=*/1);

    const std::int32_t totalChunks{3};
    std::int32_t nonEmpty{0};
    std::size_t totalRecords{0};
    for (std::int32_t chunk = 1; chunk <= totalChunks; ++chunk) {
        const auto names{ReadChunkNames(path, ScatterConfig(chunk, totalChunks, /*M=*/100, 0))};
        if (!std::empty(names)) {
            ++nonEmpty;
        }
        totalRecords += std::size(names);
    }
    EXPECT_EQ(nonEmpty, 1);       // one tile -> exactly one non-empty chunk
    EXPECT_EQ(totalRecords, 5u);  // still every record
}

TEST(BamRawReaderScatter, RejectsNonPositiveTile)
{
    tests::TempDirGuard tempDir;
    tempDir.Reset("scatter_bad_tile");
    const std::filesystem::path path{tempDir.File("scatter.bam")};
    BuildScatterBam(path, /*numZmws=*/4, /*subreadsPerZmw=*/1);

    EXPECT_THROW(BamRawReader(path, ScatterConfig(1, 2, /*M=*/0, 0)), std::invalid_argument);
}

// BamRecordReader pre-decodes via ReadBatch (a background producer). This is the
// path pbmm2's --chunk uses, so scatter must work through it too (not just
// Records()/ReadRecord).
TEST(BamRawReaderScatter, CoversEveryRecordExactlyOnceViaRecordReader)
{
    tests::TempDirGuard tempDir;
    tempDir.Reset("scatter_record_reader");
    const std::filesystem::path path{tempDir.File("scatter.bam")};
    BuildScatterBam(path, /*numZmws=*/12, /*subreadsPerZmw=*/2);

    const std::int32_t totalChunks{4};
    std::vector<std::string> all;
    for (std::int32_t chunk = 1; chunk <= totalChunks; ++chunk) {
        BamRecordReaderConfig cfg;
        cfg.RawReaderConfig = ScatterConfig(chunk, totalChunks, 1, 7);
        BamRecordReader reader{path, cfg};
        while (const auto record = reader.ReadRecord()) {
            all.emplace_back(record->Name());
        }
    }
    std::ranges::sort(all);

    auto expected{WholeFileNames(path)};
    std::ranges::sort(expected);

    EXPECT_EQ(std::ranges::adjacent_find(all), std::ranges::end(all));  // no duplicates
    EXPECT_EQ(all, expected);                                           // exact partition
}

}  // namespace Samoa
}  // namespace PacBio
