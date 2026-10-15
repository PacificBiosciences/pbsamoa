#include "Subsample.hpp"

#include "../../PathUtils.hpp"
#include "../../ZmwUtils.hpp"
#include "../CliUtils.hpp"
#include "../ParseUtils.hpp"
#include "../StatusStream.hpp"

#include <pbsamoa/PbSamoaLibraryInfo.hpp>
#include <pbsamoa/core/Bgzf.hpp>
#include <pbsamoa/core/Endian.hpp>
#include <pbsamoa/core/RawRecord.hpp>
#include <pbsamoa/core/SamHeader.hpp>
#include <pbsamoa/index/ZmwIndex.hpp>
#include <pbsamoa/io/BamRawReader.hpp>
#include <pbsamoa/io/BamWriter.hpp>

#include <pbcopper/cli2/Interface.h>
#include <pbcopper/cli2/Option.h>
#include <pbcopper/cli2/PositionalArgument.h>
#include <pbcopper/cli2/Results.h>

#include <algorithm>
#include <atomic>
#include <charconv>
#include <filesystem>
#include <format>
#include <optional>
#include <print>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <exception>

namespace PacBio {
namespace Samoa {
namespace Subsample {
namespace {

constexpr std::int32_t MAX_THREADS{64};
constexpr std::size_t BAM_TEXT_OFFSET{8};  // magic(4) + l_text(4)
constexpr std::size_t MAX_BGZF_BLOCK{65536};
// samtools keeps 24 bits of the read-name hash. The yield mode uses the same resolution.
constexpr std::uint32_t NUM_BUCKETS{std::uint32_t{1} << 24};

const CLI_v2::Option Fraction{
    R"({
    "names" : ["fraction"],
    "description" : "Fraction of templates to keep, in (0, 1]. All records with the same read name are kept or dropped together. Same records as samtools view --subsample.",
    "type" : "float",
    "default" : 0
})"};

const CLI_v2::Option Yield{
    R"({
    "names" : ["yield"],
    "description" : "Keep random ZMWs until their bases reach this yield, e.g. 60G. Suffixes K, M, G, T are powers of 1000. Reads the BAM once with a .pbi, twice without.",
    "type" : "string",
    "default" : ""
})"};

const CLI_v2::Option Seed{
    R"({
    "names" : ["seed"],
    "description" : "Sampling seed: an unsigned integer, or auto to derive it from the input header (samtools >= 1.24 default).",
    "type" : "string",
    "default" : "auto"
})"};

const CLI_v2::Option Threads{
    R"({
    "names" : ["threads"],
    "description" : "BGZF threads for reading and for writing; -1 = auto (min of hardware concurrency and 8).",
    "type" : "integer",
    "default" : -1
})"};

const CLI_v2::Option Compression{
    R"({
    "names" : ["compression"],
    "description" : "Output BGZF compression level in [1, 12].",
    "type" : "integer",
    "default" : 6
})"};

const CLI_v2::Option StatusJson{
    R"({
    "names" : ["status-json"],
    "description" : "Write a live status stream to this file: one JSON object per line (run info, phases, progress each second, final).",
    "type" : "string",
    "default" : ""
})"};

const CLI_v2::PositionalArgument Input{
    R"({
    "name" : "input",
    "description" : "Input BAM file.",
    "type" : "file"
})"};

const CLI_v2::PositionalArgument Output{
    R"({
    "name" : "output",
    "description" : "Output BAM file.",
    "type" : "file"
})"};

// The hashes below copy htslib khash.h so that the kept records match
// `samtools view --subsample`.

// __ac_FNV1a_hash_string
std::uint32_t Fnv1aHash(const std::string_view text)
{
    std::uint32_t h{2166136261};
    for (const std::uint8_t c : text) {
        h = (h ^ c) * 16777619;
    }
    return h;
}

// __ac_X31_hash_string. The char converts to khint_t as in C, so a signed char
// sign-extends on the same platforms as in samtools.
std::uint32_t X31Hash(const std::string_view text)
{
    std::uint32_t h{0};
    for (const char c : text) {
        const std::uint32_t value = c;
        h = (h << 5) - h + value;
    }
    return h;
}

// __ac_Wang_hash
std::uint32_t WangHash(std::uint32_t key)
{
    key += ~(key << 15);
    key ^= (key >> 10);
    key += (key << 3);
    key ^= (key >> 6);
    key += ~(key << 11);
    key ^= (key >> 16);
    return key;
}

/// The header text of a BAM file as stored, up to the first NUL. samtools hashes
/// these bytes for the auto seed. SamHeader::ToText() reorders lines, so it cannot
/// stand in for them.
std::string RawHeaderText(const std::filesystem::path& path)
{
    BgzfReader bgzf{path};
    std::vector<std::byte> block(MAX_BGZF_BLOCK);
    std::vector<std::byte> data;
    std::optional<std::size_t> textEnd;
    while (!textEnd || (std::size(data) < *textEnd)) {
        const std::optional<std::size_t> bytesRead{bgzf.ReadBlock(block)};
        if (!bytesRead || (*bytesRead == 0)) {
            throw std::runtime_error{
                std::format("[pbsamoa] truncated BAM header in {}", path.string())};
        }
        const std::span<const std::byte> decompressed{std::data(block), *bytesRead};
        data.insert(std::end(data), std::begin(decompressed), std::end(decompressed));
        if (!textEnd && (std::size(data) >= BAM_TEXT_OFFSET)) {
            constexpr std::string_view BAM_MAGIC{"BAM\1", 4};
            const std::string_view magic{reinterpret_cast<const char*>(std::data(data)), 4};
            if (magic != BAM_MAGIC) {
                throw std::runtime_error{
                    std::format("[pbsamoa] not a BAM file: {}", path.string())};
            }
            textEnd = BAM_TEXT_OFFSET + ReadU32LE(std::data(data) + 4);
        }
    }
    const std::string_view text{reinterpret_cast<const char*>(std::data(data)) + BAM_TEXT_OFFSET,
                                *textEnd - BAM_TEXT_OFFSET};
    return std::string{text.substr(0, text.find('\0'))};
}

/// The seed before samtools scrambles it. This is the value in the @CO line.
std::uint32_t ResolveSeed(const std::string_view seedText, const std::filesystem::path& input)
{
    if (seedText == "auto") {
        return Fnv1aHash(RawHeaderText(input));
    }
    return Tools::ParseIntegerOrThrow<std::uint32_t>(seedText, "--seed");
}

/// samtools scrambles a non-zero seed with the C library rand(). glibc and macOS
/// libc have different rand() functions, so the kept records match samtools on the
/// same platform only.
std::uint32_t ScrambleSeed(const std::uint32_t seed)
{
    if (seed == 0) {
        return 0;
    }
    std::srand(seed);
    return std::rand();
}

/// The 24-bit read-name bucket of samtools. samtools keeps a record when
/// bucket / 2^24 < fraction.
std::uint32_t NameBucket(const std::string_view name, const std::uint32_t scrambledSeed)
{
    return WangHash(X31Hash(name) ^ scrambledSeed) & (NUM_BUCKETS - 1);
}

/// bucket / 2^24 < fraction is the same as bucket < ceil(fraction * 2^24), because the
/// bucket is an integer and the scale by 2^24 is exact.
std::uint32_t CutoffForFraction(const double fraction)
{
    const std::uint32_t cutoff = std::ceil(fraction * NUM_BUCKETS);
    return cutoff;
}

/// The 24-bit bucket of a ZMW for the yield mode: the top bits of a SplitMix64 mix of
/// (rgId, zmw) and the seed. A .pbi has these two values but no read names. The mix
/// does not use rand(), so the selection is the same on every platform.
std::uint32_t ZmwBucket(const ZmwIdentity id, const std::uint32_t seed)
{
    const std::uint32_t rgBits = id.rgId;
    const std::uint32_t zmwBits = id.zmw;
    std::uint64_t x{((std::uint64_t{rgBits} << 32) | zmwBits) ^
                    (std::uint64_t{seed} * 0x9E3779B97F4A7C15)};
    x ^= x >> 30;
    x *= 0xBF58476D1CE4E5B9;
    x ^= x >> 27;
    x *= 0x94D049BB133111EB;
    x ^= x >> 31;
    return x >> 40;
}

/// The (rgId, zmw) identity of a record for the yield mode. The name must be a PacBio
/// name movie/zmw/..., because a name without a ZMW would put every record in one
/// bucket.
ZmwIdentity YieldIdentity(const RawRecordView& view)
{
    const std::string_view name{view.Name()};
    const std::size_t firstSlash{name.find('/')};
    const std::size_t zmwStart{firstSlash == std::string_view::npos ? std::size(name)
                                                                    : firstSlash + 1};
    const std::string_view zmwField{name.substr(zmwStart, name.find('/', zmwStart) - zmwStart)};
    if (!Tools::ParseInteger<std::int32_t>(zmwField, "ZMW")) {
        throw std::runtime_error{
            std::format("--yield needs PacBio read names movie/zmw/..., found: {}", name)};
    }
    return detail::ParseRecordIdentity(view.RawData());
}

/// Parse a yield in bases with an optional decimal suffix: 60G is 60 * 10^9 bases.
std::int64_t ParseYield(const std::string_view text)
{
    // Suffix pair i (lower and upper case) scales by 1000^(i + 1).
    constexpr std::string_view SUFFIXES{"kKmMgGtT"};
    std::string_view number{text};
    double scale{1.0};
    if (!std::empty(number)) {
        if (const std::size_t pos{SUFFIXES.find(number.back())}; pos != std::string_view::npos) {
            scale = std::pow(1000.0, (pos / 2) + 1);
            number.remove_suffix(1);
        }
    }
    double value{0.0};
    const char* const end{std::data(number) + std::size(number)};
    const auto [ptr, ec]{std::from_chars(std::data(number), end, value)};
    const double bases{value * scale};
    // 9e18 stays below the std::int64_t maximum of about 9.22e18.
    if ((ec != std::errc{}) || (ptr != end) || !(bases >= 1.0) || (bases > 9e18)) {
        throw std::runtime_error{std::format(
            "invalid --yield: {} (give bases with an optional K, M, G or T suffix, e.g. 60G)",
            text)};
    }
    return std::llround(bases);
}

/// Bases per ZMW bucket, and the totals.
struct YieldHistogram
{
    std::vector<std::int64_t> BasesPerBucket = std::vector<std::int64_t>(NUM_BUCKETS);
    std::int64_t InputBases{0};
    std::int64_t InputRecords{0};

    void Add(const std::uint32_t bucket, const std::int64_t bases)
    {
        BasesPerBucket[bucket] += bases;
        InputBases += bases;
        ++InputRecords;
    }
};

/// One pass over the .pbi only. For unaligned HiFi, qEnd - qStart is the read length.
YieldHistogram HistogramFromPbi(const ZmwIndex& index, const std::uint32_t seed)
{
    YieldHistogram histogram;
    const std::span<const std::int32_t> lengths{index.QueryLengths()};
    for (std::size_t row{0}; row < std::size(lengths); ++row) {
        histogram.Add(ZmwBucket(index.Identity(row), seed), lengths[row]);
    }
    return histogram;
}

/// The extra pass over the BAM when there is no .pbi. Adds the records of each batch to
/// \p done for the status stream.
YieldHistogram HistogramFromBam(const std::filesystem::path& input, const std::uint32_t seed,
                                const std::size_t numThreads, std::atomic<std::int64_t>& done)
{
    YieldHistogram histogram;
    BamRawReader reader{input, BamRawReaderConfig{.BgzfWorkers = numThreads}};
    while (const std::optional<RawRecordBatch> batch{reader.ReadBatch()}) {
        for (std::size_t i{0}; i < batch->RecordCount(); ++i) {
            const RawRecordView view{batch->View(i)};
            histogram.Add(ZmwBucket(YieldIdentity(view), seed), view.SeqLength());
        }
        done.fetch_add(batch->RecordCount(), std::memory_order::relaxed);
    }
    return histogram;
}

/// The smallest cutoff whose buckets [0, cutoff) hold at least targetBases. The kept
/// yield exceeds the target by less than the bases of the last bucket.
std::uint32_t CutoffForYield(const YieldHistogram& histogram, const std::int64_t targetBases)
{
    std::int64_t bases{0};
    for (std::uint32_t bucket{0}; bucket < NUM_BUCKETS; ++bucket) {
        bases += histogram.BasesPerBucket[bucket];
        if (bases >= targetBases) {
            return bucket + 1;
        }
    }
    return NUM_BUCKETS;
}

/// Throws when record \p row of the BAM is not row \p row of its .pbi. A stale index
/// would give a wrong cutoff.
void CheckPbiRow(const ZmwIndex& index, const std::size_t row, const ZmwIdentity id,
                 const std::uint32_t seqLength)
{
    if ((row >= index.NumRecords()) || (index.Identity(row) != id) ||
        std::cmp_not_equal(index.QueryLengths()[row], seqLength)) {
        throw std::runtime_error{std::format(
            "[pbsamoa] the .pbi does not match BAM record {}; rebuild it with pbindex", row)};
    }
}

struct Counts
{
    std::int64_t Records{0};
    std::int64_t Bases{0};
    std::int64_t KeptRecords{0};
    std::int64_t KeptBases{0};
};

/// The counts as the metrics of the status stream.
JSON::Json CountsJson(const Counts& counts)
{
    return {{"records", counts.Records},
            {"bases", counts.Bases},
            {"kept_records", counts.KeptRecords},
            {"kept_bases", counts.KeptBases}};
}

/// Copy each record whose bucket is below \p cutoff. Adds the records of each batch to
/// \p done for the status stream.
template <typename BUCKET_OF>
Counts WriteSubsample(BamRawReader& reader, BamWriter& writer, const std::uint32_t cutoff,
                      const BUCKET_OF& bucketOf, std::atomic<std::int64_t>& done)
{
    Counts counts;
    while (const std::optional<RawRecordBatch> batch{reader.ReadBatch()}) {
        for (std::size_t i{0}; i < batch->RecordCount(); ++i) {
            const RawRecordView view{batch->View(i)};
            const std::int64_t bases{view.SeqLength()};
            ++counts.Records;
            counts.Bases += bases;
            if (bucketOf(view) < cutoff) {
                writer.Write(view);
                ++counts.KeptRecords;
                counts.KeptBases += bases;
            }
        }
        done.fetch_add(batch->RecordCount(), std::memory_order::relaxed);
    }
    return counts;
}

struct Settings
{
    std::filesystem::path Input;
    std::filesystem::path Output;
    double Fraction{0.0};
    std::optional<std::int64_t> TargetBases;  // set in the yield mode
    std::uint32_t Seed{0};
    std::size_t NumThreads{0};
    std::int32_t Compression{6};
    std::filesystem::path StatusJson;
};

Settings ParseSettings(const CLI_v2::Results& results)
{
    const std::vector<std::string>& positional{results.PositionalArguments()};
    if (std::size(positional) != 2) {
        throw std::runtime_error{
            "subsample requires exactly two arguments: <input.bam> <output.bam>"};
    }
    Settings settings{.Input = positional[0], .Output = positional[1]};

    settings.Fraction = results[Fraction];
    const std::string yieldText = results[Yield];
    if (!std::empty(yieldText)) {
        if (settings.Fraction != 0.0) {
            throw std::runtime_error{"give --fraction or --yield, not both"};
        }
        settings.TargetBases = ParseYield(yieldText);
    } else if (!(settings.Fraction > 0.0) || (settings.Fraction > 1.0)) {
        // samtools ignores --subsample 0 and keeps every record; a fraction of 0 is an
        // error here.
        throw std::runtime_error{
            std::format("--fraction must be in (0, 1]: {}", settings.Fraction)};
    }

    const std::string seedText = results[Seed];
    settings.Seed = ResolveSeed(seedText, settings.Input);
    const std::int32_t requestedThreads{results[Threads]};
    settings.NumThreads = Tools::ResolveNumWorkers(requestedThreads, MAX_THREADS);
    const std::int32_t compression{results[Compression]};
    settings.Compression = Tools::CheckCompressionLevel(compression);
    const std::string statusJson = results[StatusJson];
    settings.StatusJson = statusJson;
    return settings;
}

struct YieldPlan
{
    std::optional<ZmwIndex> Pbi;
    std::uint32_t Cutoff{0};
    // the records of the input: the total of the write phase
    std::int64_t InputRecords{0};
};

/// The .pbi of the input, if it exists.
std::filesystem::path PbiPath(const Settings& settings)
{
    return SidecarPath(settings.Input, ".pbi");
}

/// The yield mode reads the lengths from the .pbi when there is one, so that it reads
/// the BAM only once (phase `select`). Without a .pbi, an extra pass over the BAM gets
/// them (phase `count`).
YieldPlan PlanYield(const Settings& settings, Tools::StatusStream& status)
{
    YieldPlan plan;
    const bool hasPbi{std::filesystem::exists(PbiPath(settings))};
    std::atomic<std::int64_t> counted{0};
    Tools::StatusStream::Phase phase{&status, hasPbi ? "select" : "count",
                                     hasPbi ? nullptr : &counted, "reads"};
    if (hasPbi) {
        plan.Pbi = ZmwIndex::FromPbi(PbiPath(settings));
    }
    const YieldHistogram histogram{
        plan.Pbi ? HistogramFromPbi(*plan.Pbi, settings.Seed)
                 : HistogramFromBam(settings.Input, settings.Seed, settings.NumThreads, counted)};
    plan.Cutoff = CutoffForYield(histogram, *settings.TargetBases);
    plan.InputRecords = histogram.InputRecords;
    phase.End({{"records", histogram.InputRecords},
               {"bases", histogram.InputBases},
               {"target_bases", *settings.TargetBases}});
    if (histogram.InputBases < *settings.TargetBases) {
        std::println(stderr,
                     "pbsamoa subsample: WARNING: the input has {} bases, less than the target "
                     "yield {}; keeping all records",
                     histogram.InputBases, *settings.TargetBases);
    }
    return plan;
}

Counts WriteYieldSubsample(BamRawReader& reader, BamWriter& writer, const YieldPlan& plan,
                           const std::uint32_t seed, std::atomic<std::int64_t>& done)
{
    std::size_t row{0};
    const Counts counts{WriteSubsample(
        reader, writer, plan.Cutoff,
        [&](const RawRecordView& view) {
            const ZmwIdentity id{YieldIdentity(view)};
            if (plan.Pbi) {
                CheckPbiRow(*plan.Pbi, row, id, view.SeqLength());
                ++row;
            }
            return ZmwBucket(id, seed);
        },
        done)};
    if (plan.Pbi && (row != plan.Pbi->NumRecords())) {
        throw std::runtime_error{
            std::format("[pbsamoa] the .pbi has {} records, the BAM {}; rebuild it with pbindex",
                        plan.Pbi->NumRecords(), row)};
    }
    return counts;
}

Counts WriteFractionSubsample(BamRawReader& reader, BamWriter& writer, const double fraction,
                              const std::uint32_t seed, std::atomic<std::int64_t>& done)
{
    const std::uint32_t scrambledSeed{ScrambleSeed(seed)};
    return WriteSubsample(
        reader, writer, CutoffForFraction(fraction),
        [&](const RawRecordView& view) { return NameBucket(view.Name(), scrambledSeed); }, done);
}

/// Open the output with \p comment and a @PG record, and fill it with \p write (phase
/// `write`, in reads up to \p total). On any error, remove the output, so that no partial
/// or wrong BAM stays behind.
template <typename WRITE>
Counts WriteOutput(const Settings& settings, const std::string& commandLine,
                   const std::string& comment, Tools::StatusStream& status,
                   const std::optional<std::int64_t> total, const WRITE& write)
{
    BamRawReader reader{settings.Input, BamRawReaderConfig{.BgzfWorkers = settings.NumThreads}};
    SamHeader header{reader.Header()};
    header.AddComment(comment);
    AppendProgramRecord(header, "pbsamoa.subsample", commandLine);
    try {
        BamWriter writer{
            settings.Output, header,
            BamWriterConfig{
                .BgzfConfig =
                    {
                        .CompressionLevel = settings.Compression,
                        .BgzfWorkers = std::ranges::max(settings.NumThreads, std::size_t{1}),
                    },
            }};
        std::atomic<std::int64_t> done{0};
        Tools::StatusStream::Phase phase{&status, "write", &done, "reads", total};
        const Counts counts{write(reader, writer, done)};
        writer.Close();
        phase.End(CountsJson(counts));
        return counts;
    } catch (...) {
        std::error_code ignored;
        std::filesystem::remove(settings.Output, ignored);
        throw;
    }
}

void PrintCounts(const Counts& counts, const std::string_view details)
{
    std::println(stderr, "pbsamoa subsample: kept {} of {} records, {} of {} bases ({})",
                 counts.KeptRecords, counts.Records, counts.KeptBases, counts.Bases, details);
}

}  // namespace

CLI_v2::Interface CreateInterface()
{
    CLI_v2::Interface interface{
        "pbsamoa subsample",
        "Keep a fraction of the templates in a BAM file (same records as samtools view "
        "--subsample), or random ZMWs up to a target yield.",
        LibraryFormattedVersion()};
    // `subsample` keeps its own --threads (-1 = auto, capped at min(hw, 8)); the built-in
    // --num-threads resolves 0 to the raw hardware count, which would drop that cap.
    interface.DisableNumThreadsOption();
    interface.AddOptions({Fraction, Yield, Seed, Threads, Compression, StatusJson});
    interface.AddPositionalArguments({Input, Output});
    return interface;
}

int Runner(const CLI_v2::Results& results)
{
    const Settings settings{ParseSettings(results)};
    const std::string commandLine{"pbsamoa " + results.InputCommandLine()};
    // outside the try: the catch writes the final line of a failure
    Tools::StatusStream status;
    std::vector<std::string> phases{"write"};
    if (settings.TargetBases) {
        phases.insert(std::begin(phases),
                      std::filesystem::exists(PbiPath(settings)) ? "select" : "count");
    }
    status.Open(settings.StatusJson,
                Tools::MakeRunInfo("pbsamoa", "subsample", LibraryFormattedVersion(), commandLine,
                                   settings.NumThreads, phases));
    try {
        Counts counts;
        JSON::Json metrics;
        if (settings.TargetBases) {
            const YieldPlan plan{PlanYield(settings, status)};
            counts = WriteOutput(
                settings, commandLine,
                std::format("Sub-sampled yield={} seed={}", *settings.TargetBases, settings.Seed),
                status, plan.InputRecords,
                [&](BamRawReader& reader, BamWriter& writer, std::atomic<std::int64_t>& done) {
                    return WriteYieldSubsample(reader, writer, plan, settings.Seed, done);
                });
            PrintCounts(counts,
                        std::format("target {}, seed {}, {}", *settings.TargetBases, settings.Seed,
                                    plan.Pbi ? "1 BAM pass with .pbi" : "2 BAM passes, no .pbi"));
            metrics = CountsJson(counts);
            metrics["target_bases"] = *settings.TargetBases;
        } else {
            // the .pbi gives the total of the write phase cheaply, if there is one
            std::optional<std::int64_t> total;
            if (std::filesystem::exists(PbiPath(settings))) {
                const std::int64_t numRecords = ZmwIndex::FromPbi(PbiPath(settings)).NumRecords();
                total = numRecords;
            }
            counts = WriteOutput(
                settings, commandLine,
                std::format("Sub-sampled fraction={:g} seed={}", settings.Fraction, settings.Seed),
                status, total,
                [&](BamRawReader& reader, BamWriter& writer, std::atomic<std::int64_t>& done) {
                    return WriteFractionSubsample(reader, writer, settings.Fraction, settings.Seed,
                                                  done);
                });
            PrintCounts(counts,
                        std::format("fraction {:g}, seed {}", settings.Fraction, settings.Seed));
            metrics = CountsJson(counts);
            metrics["fraction"] = settings.Fraction;
        }
        metrics["seed"] = settings.Seed;
        status.Final(
            "FINISHED", {}, std::move(metrics),
            JSON::Json::array({{{"path", std::filesystem::absolute(settings.Output).string()},
                                {"kind", "bam"}}}));
    } catch (const std::exception& e) {
        status.Final("FAILED", e.what());
        throw;
    }
    return EXIT_SUCCESS;
}

}  // namespace Subsample
}  // namespace Samoa
}  // namespace PacBio
