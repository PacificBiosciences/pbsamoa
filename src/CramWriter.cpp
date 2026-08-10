#include <pbsamoa/io/CramWriter.hpp>

#include <pbsamoa/cram/CramCodec.hpp>
#include <pbsamoa/cram/CramCompression.hpp>
#include <pbsamoa/cram/CramStructs.hpp>
#include <pbsamoa/index/CraiIndex.hpp>

#include "BinaryUtils.hpp"
#include "CramInternal.hpp"
#include "CramStructs.hpp"
#include "WriterUtils.hpp"

#include <pbcopper/parallel/ThreadPool.h>

#include <htscodecs/fqzcomp_qual.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <format>
#include <fstream>
#include <future>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace PacBio {
namespace Samoa {

// Defined in CramStructs.cpp; fills Landmarks, Length, and NumBlocks from the container.
CramContainerHeader ComputeContainerLandmarks(const CramContainer& container,
                                              std::int64_t compressionHeaderBlockSize);

namespace CramWriterInternal {

struct EncodedTagPayload
{
    std::int32_t ContentId{};
    char Type{};
    std::vector<std::byte> Data;
};

struct RecordTagData
{
    std::int32_t ReadGroupIndex{-1};
    std::int32_t TagListIndex{0};
    std::vector<TagTriple> TagSet;
    std::vector<EncodedTagPayload> Payloads;
};

struct RefSpan
{
    std::int32_t SequenceId{};
    bool HasMappedSpan{false};
    std::int64_t MinPos{};
    std::int64_t MaxEnd{};
};

RecordTagData BuildRecordTagData(const BamRecord& record)
{
    RecordTagData result;

    for (const auto& [key, value] : record.Tags().Entries()) {
        if (key == RG_TAG) {
            continue;
        }

        auto encoded = EncodeTagValueToBamPayload(value);
        const auto type = encoded.Type;
        result.TagSet.push_back(TagTriple{key.First(), key.Second(), type});
        result.Payloads.push_back(EncodedTagPayload{
            .ContentId = TagContentId(key, type),
            .Type = type,
            .Data = std::move(encoded.Payload),
        });
    }

    std::ranges::sort(result.TagSet);

    return result;
}

std::vector<std::byte> BuildTagIdsDictionary(const std::vector<std::vector<TagTriple>>& tagSets)
{
    std::vector<std::byte> dict;
    for (const auto& tagSet : tagSets) {
        for (const auto& triple : tagSet) {
            dict.push_back(static_cast<std::byte>(triple.Tag1));
            dict.push_back(static_cast<std::byte>(triple.Tag2));
            dict.push_back(static_cast<std::byte>(triple.Type));
        }
        dict.push_back(std::byte{0});
    }
    return dict;
}

struct DataSeriesInfo
{
    CramDataSeries Series;
    std::int32_t BlockId;
};

constexpr std::array<DataSeriesInfo, 28> DATA_SERIES = {{
    {CramDataSeries::BF, 1},  {CramDataSeries::CF, 2},  {CramDataSeries::RI, 3},
    {CramDataSeries::RL, 4},  {CramDataSeries::AP, 5},  {CramDataSeries::RG, 6},
    {CramDataSeries::RN, 7},  {CramDataSeries::MF, 8},  {CramDataSeries::NS, 9},
    {CramDataSeries::NP, 10}, {CramDataSeries::TS, 11}, {CramDataSeries::NF, 12},
    {CramDataSeries::TL, 13}, {CramDataSeries::FN, 14}, {CramDataSeries::FC, 15},
    {CramDataSeries::FP, 16}, {CramDataSeries::MQ, 17}, {CramDataSeries::BA, 18},
    {CramDataSeries::QS, 19}, {CramDataSeries::BS, 20}, {CramDataSeries::IN, 21},
    {CramDataSeries::DL, 22}, {CramDataSeries::SC, 23}, {CramDataSeries::RS, 24},
    {CramDataSeries::PD, 25}, {CramDataSeries::HC, 26}, {CramDataSeries::BB, 27},
    {CramDataSeries::QQ, 28},
}};

constexpr std::int32_t BlockIdFor(CramDataSeries ds)
{
    for (const auto& info : DATA_SERIES) {
        if (info.Series == ds) {
            return info.BlockId;
        }
    }
    return -1;
}

constexpr std::int32_t BF_BLOCK_ID = BlockIdFor(CramDataSeries::BF);
constexpr std::int32_t CF_BLOCK_ID = BlockIdFor(CramDataSeries::CF);
constexpr std::int32_t RI_BLOCK_ID = BlockIdFor(CramDataSeries::RI);
constexpr std::int32_t RL_BLOCK_ID = BlockIdFor(CramDataSeries::RL);
constexpr std::int32_t AP_BLOCK_ID = BlockIdFor(CramDataSeries::AP);
constexpr std::int32_t RG_BLOCK_ID = BlockIdFor(CramDataSeries::RG);
constexpr std::int32_t RN_BLOCK_ID = BlockIdFor(CramDataSeries::RN);
constexpr std::int32_t MF_BLOCK_ID = BlockIdFor(CramDataSeries::MF);
constexpr std::int32_t NS_BLOCK_ID = BlockIdFor(CramDataSeries::NS);
constexpr std::int32_t NP_BLOCK_ID = BlockIdFor(CramDataSeries::NP);
constexpr std::int32_t TS_BLOCK_ID = BlockIdFor(CramDataSeries::TS);
constexpr std::int32_t TL_BLOCK_ID = BlockIdFor(CramDataSeries::TL);
constexpr std::int32_t FN_BLOCK_ID = BlockIdFor(CramDataSeries::FN);
constexpr std::int32_t FC_BLOCK_ID = BlockIdFor(CramDataSeries::FC);
constexpr std::int32_t FP_BLOCK_ID = BlockIdFor(CramDataSeries::FP);
constexpr std::int32_t MQ_BLOCK_ID = BlockIdFor(CramDataSeries::MQ);
constexpr std::int32_t BA_BLOCK_ID = BlockIdFor(CramDataSeries::BA);
constexpr std::int32_t QS_BLOCK_ID = BlockIdFor(CramDataSeries::QS);
constexpr std::int32_t BB_BLOCK_ID = BlockIdFor(CramDataSeries::BB);
constexpr std::int32_t QQ_BLOCK_ID = BlockIdFor(CramDataSeries::QQ);
constexpr std::int32_t IN_BLOCK_ID = BlockIdFor(CramDataSeries::IN);
constexpr std::int32_t DL_BLOCK_ID = BlockIdFor(CramDataSeries::DL);
constexpr std::int32_t SC_BLOCK_ID = BlockIdFor(CramDataSeries::SC);
constexpr std::int32_t RS_BLOCK_ID = BlockIdFor(CramDataSeries::RS);
constexpr std::int32_t PD_BLOCK_ID = BlockIdFor(CramDataSeries::PD);
constexpr std::int32_t HC_BLOCK_ID = BlockIdFor(CramDataSeries::HC);

RefSpan& LookupOrCreateRefSpan(std::vector<RefSpan>& spans,
                               std::unordered_map<std::int32_t, std::size_t>& spanLookup,
                               std::int32_t sequenceId)
{
    if (const auto it = spanLookup.find(sequenceId); it != std::end(spanLookup)) {
        return spans[it->second];
    }

    const std::size_t idx = std::size(spans);
    spans.push_back(RefSpan{.SequenceId = sequenceId});
    spanLookup.emplace(sequenceId, idx);
    return spans.back();
}

CramEncodingDescriptor MakeExternalDescriptor(std::int32_t blockId)
{
    CramEncodingDescriptor desc;
    desc.CodecId = CramCodecId::EXTERNAL;
    WriteItf8(desc.Parameters, blockId);
    return desc;
}

void AppendNestedEncodingDescriptor(std::vector<std::byte>& parameters,
                                    const CramEncodingDescriptor& desc)
{
    WriteItf8(parameters, std::to_underlying(desc.CodecId));
    WriteItf8(parameters, static_cast<std::int32_t>(std::size(desc.Parameters)));
    parameters.insert(std::end(parameters), std::begin(desc.Parameters), std::end(desc.Parameters));
}

CramEncodingDescriptor MakeByteArrayLenExternalDescriptor(std::int32_t blockId)
{
    const auto lenDesc = MakeExternalDescriptor(blockId);
    const auto dataDesc = MakeExternalDescriptor(blockId);

    CramEncodingDescriptor desc;
    desc.CodecId = CramCodecId::BYTE_ARRAY_LEN;
    AppendNestedEncodingDescriptor(desc.Parameters, lenDesc);
    AppendNestedEncodingDescriptor(desc.Parameters, dataDesc);
    return desc;
}

CramEncodingDescriptor MakeByteArrayStopExternalDescriptor(std::int32_t blockId, std::byte stopByte)
{
    CramEncodingDescriptor desc;
    desc.CodecId = CramCodecId::BYTE_ARRAY_STOP;
    desc.Parameters.push_back(stopByte);
    WriteItf8(desc.Parameters, blockId);
    return desc;
}

CramEncodingDescriptor MakeTagEncodingDescriptor(char type, std::int32_t blockId)
{
    if (type == 'Z' || type == 'H') {
        return MakeByteArrayStopExternalDescriptor(blockId, std::byte{'\t'});
    }
    return MakeByteArrayLenExternalDescriptor(blockId);
}

CramEncodingDescriptor MakeDataSeriesEncodingDescriptor(CramDataSeries series, std::int32_t blockId)
{
    switch (series) {
        case CramDataSeries::RN:
        case CramDataSeries::IN:
        case CramDataSeries::SC:
            return MakeByteArrayStopExternalDescriptor(blockId, std::byte{0});
        case CramDataSeries::BB:
            return MakeByteArrayLenExternalDescriptor(blockId);
        default:
            return MakeExternalDescriptor(blockId);
    }
}

std::string DataSeriesCode(CramDataSeries ds)
{
    const auto value = std::to_underlying(ds);
    std::string code(2, '\0');
    code[0] = static_cast<char>((value >> 8) & 0xFF);
    code[1] = static_cast<char>(value & 0xFF);
    return code;
}

std::int64_t ReferenceEndOrNextBase(const BamRecord& record)
{
    const std::int64_t pos{record.Pos()};
    return NonEmptyAlignmentEnd(pos, static_cast<std::int64_t>(record.ReferenceEnd()));
}

void UpdateRefSpan(RefSpan& span, const BamRecord& record)
{
    const std::int64_t start = record.Pos();
    const std::int64_t end = ReferenceEndOrNextBase(record);
    if (!span.HasMappedSpan) {
        span.HasMappedSpan = true;
        span.MinPos = start;
        span.MaxEnd = end;
        return;
    }

    span.MinPos = std::min(span.MinPos, start);
    span.MaxEnd = std::max(span.MaxEnd, end);
}

std::int32_t RecordAlignmentEnd(const BamRecord& record)
{
    if (record.IsMapped()) {
        return static_cast<std::int32_t>(ReferenceEndOrNextBase(record));
    }
    return record.Pos() + static_cast<std::int32_t>(std::size(record.Sequence()));
}

bool IsV31BlockMethod(const CramBlockMethod method)
{
    switch (method) {
        case CramBlockMethod::RANS4X16:
        case CramBlockMethod::ADAPTIVE_ARITH:
        case CramBlockMethod::FQZCOMP:
        case CramBlockMethod::NAME_TOKENISER:
            return true;
        default:
            return false;
    }
}

bool DataSeriesSupportsFqzcomp(const CramDataSeries dataSeries)
{
    return (dataSeries == CramDataSeries::QS) || (dataSeries == CramDataSeries::QQ);
}

libdeflate_compressor* ThreadLocalGzipCompressor(const std::optional<int> compressionLevel)
{
    const int level = compressionLevel.value_or(DEFAULT_GZIP_COMPRESSION_LEVEL);
    if (level < 0 || level > 12) {
        throw std::runtime_error(
            std::format("CramWriter: invalid gzip compression level {}", level));
    }

    struct ThreadLocalCompressorCache
    {
        int Level{-1};
        LibdeflateCompressorPtr Compressor{};
    };

    static thread_local ThreadLocalCompressorCache cache;

    if (!cache.Compressor || cache.Level != level) {
        cache.Compressor.reset(libdeflate_alloc_compressor(level));
        cache.Level = level;
    }

    if (!cache.Compressor) {
        throw std::runtime_error("CramWriter: failed to allocate gzip compressor");
    }
    return cache.Compressor.get();
}

CramBlockMethod ResolveBlockMethod(CramBlockContentType contentType, std::int32_t contentId,
                                   const std::unordered_map<std::int32_t, CramBlockMethod>& methods,
                                   CramBlockMethod fallbackMethod)
{
    if (contentType != CramBlockContentType::EXTERNAL_DATA) {
        return fallbackMethod;
    }
    if (const auto it = methods.find(contentId); it != std::end(methods)) {
        return it->second;
    }
    return fallbackMethod;
}

libdeflate_compressor* SelectGzipCompressor(bool useThreadLocalGzip,
                                            LibdeflateCompressorPtr& gzipCompressor,
                                            std::optional<int> compressionLevel)
{
    if (useThreadLocalGzip) {
        return ThreadLocalGzipCompressor(compressionLevel);
    }
    return gzipCompressor.get();
}

void CompressFqzDataBlock(CramBlock& block, std::uint64_t totalQualityBytes,
                          std::span<const std::uint32_t> qualityRecordLengths,
                          std::span<const std::uint32_t> qualityRecordFlags)
{
    if (totalQualityBytes != static_cast<std::uint64_t>(block.RawSize)) {
        throw std::runtime_error(
            std::format("CramWriter: fqz metadata size mismatch for content ID {} "
                        "(metadata={}, raw={})",
                        block.ContentId, totalQualityBytes, block.RawSize));
    }

    block.Data = CramFqzcompCompress(block.Data, qualityRecordLengths, qualityRecordFlags);
    block.Method = CramBlockMethod::FQZCOMP;
    if (std::size(block.Data) > std::numeric_limits<std::int32_t>::max()) {
        throw std::runtime_error("CramWriter: fqz compressed block exceeds CRAM size limit");
    }
    block.CompressedSize = static_cast<std::int32_t>(std::size(block.Data));
}

void CompressDataBlock(CramBlock& block, CramBlockMethod method, bool useThreadLocalCompressor,
                       LibdeflateCompressorPtr& gzipCompressor, std::optional<int> compressionLevel,
                       std::uint64_t totalQualityBytes,
                       std::span<const std::uint32_t> qualityRecordLengths,
                       std::span<const std::uint32_t> qualityRecordFlags)
{
    if (method == CramBlockMethod::FQZCOMP &&
        (block.ContentId == QS_BLOCK_ID || block.ContentId == QQ_BLOCK_ID)) {
        CompressFqzDataBlock(block, totalQualityBytes, qualityRecordLengths, qualityRecordFlags);
        return;
    }

    auto* compressor =
        SelectGzipCompressor(useThreadLocalCompressor, gzipCompressor, compressionLevel);
    CompressCramBlock(block, method, compressor, compressionLevel);
}

std::int32_t QueryCopyLength(std::int32_t readLength, std::int32_t requested,
                             std::int32_t seqOffset)
{
    if (requested <= 0 || seqOffset < 0 || seqOffset >= readLength) {
        return 0;
    }
    return std::ranges::min(requested, readLength - seqOffset);
}

std::span<const std::byte> QuerySpan(std::string_view seq, std::int32_t seqOffset,
                                     std::int32_t requested)
{
    const auto copyLen =
        QueryCopyLength(static_cast<std::int32_t>(std::size(seq)), requested, seqOffset);
    if (copyLen <= 0) {
        return {};
    }
    const auto* begin = reinterpret_cast<const std::byte*>(seq.data() + seqOffset);
    return {begin, static_cast<std::size_t>(copyLen)};
}

void WriteByteArrayLen(CramExternalBlockStore& extStore, std::int32_t blockId,
                       std::span<const std::byte> data)
{
    extStore.WriteItf8(blockId, static_cast<std::int32_t>(std::size(data)));
    extStore.WriteBytes(blockId, data);
}

void WriteByteArrayStop(CramExternalBlockStore& extStore, std::int32_t blockId,
                        std::span<const std::byte> data, std::byte stop)
{
    extStore.WriteBytes(blockId, data);
    extStore.WriteByte(blockId, stop);
}

void WriteEncodedTagPayload(CramExternalBlockStore& extStore, const EncodedTagPayload& payload)
{
    if (payload.Type == 'Z' || payload.Type == 'H') {
        WriteByteArrayStop(extStore, payload.ContentId, payload.Data, std::byte{'\t'});
        return;
    }
    WriteByteArrayLen(extStore, payload.ContentId, payload.Data);
}

void WriteQualities(CramExternalBlockStore& extStore, std::span<const std::uint8_t> qualities,
                    std::int32_t expectedLength)
{
    if (static_cast<std::int32_t>(std::size(qualities)) == expectedLength) {
        const auto* qPtr = reinterpret_cast<const std::byte*>(qualities.data());
        extStore.WriteBytes(QS_BLOCK_ID, {qPtr, std::size(qualities)});
        return;
    }
    for (std::int32_t qi = 0; qi < expectedLength; ++qi) {
        std::uint8_t q{0xFF};
        if (qi < static_cast<std::int32_t>(std::size(qualities))) {
            q = qualities[qi];
        }
        extStore.WriteByte(QS_BLOCK_ID, static_cast<std::byte>(q));
    }
}

void WriteFeatureHeader(CramExternalBlockStore& extStore, char code, std::int32_t featurePos,
                        std::int32_t& previousFeaturePos)
{
    extStore.WriteByte(FC_BLOCK_ID, static_cast<std::byte>(static_cast<std::uint8_t>(code)));
    extStore.WriteItf8(FP_BLOCK_ID, featurePos - previousFeaturePos);
    previousFeaturePos = featurePos;
}

}  // namespace CramWriterInternal

using namespace CramWriterInternal;

// ---------------------------------------------------------------------------
// Impl
// ---------------------------------------------------------------------------

struct CramWriter::Impl
{
    std::filesystem::path path;       // final output path requested by caller
    std::filesystem::path writePath;  // actual on-disk output (temp path when UseTempFile=true)
    std::ofstream file;
    SamHeader header;
    CramWriterConfig config;
    bool closed{false};

    std::vector<BamRecord> pendingRecords;
    std::vector<std::vector<BamRecord>> pendingSlices;
    std::int64_t globalRecordCounter{0};
    std::int64_t nextContainerOffset{0};
    std::unordered_map<std::string, std::int32_t> readGroupIndex;
    std::vector<CraiEntry> craiEntries;
    std::future<void> pendingWrite;
    LibdeflateCompressorPtr gzipCompressor;
    std::shared_ptr<Parallel::ThreadPool<>> compressionPool;

    CramBlockMethod DefaultDataSeriesMethod(CramDataSeries ds) const
    {
        switch (config.BlockCompressionMethod) {
            case CramBlockMethod::NAME_TOKENISER:
                if (ds == CramDataSeries::RN) {
                    return CramBlockMethod::NAME_TOKENISER;
                }
                return CramBlockMethod::GZIP;
            case CramBlockMethod::FQZCOMP:
                if (ds == CramDataSeries::QS || ds == CramDataSeries::QQ) {
                    return CramBlockMethod::FQZCOMP;
                }
                return CramBlockMethod::GZIP;
            default:
                return config.BlockCompressionMethod;
        }
    }

    CramBlockMethod ResolveDataSeriesMethod(CramDataSeries ds) const
    {
        CramBlockMethod method = DefaultDataSeriesMethod(ds);
        if (const auto it = config.DataSeriesCompressionMethods.find(ds);
            it != std::end(config.DataSeriesCompressionMethods)) {
            method = it->second;
        }

        if ((method == CramBlockMethod::NAME_TOKENISER) && (ds != CramDataSeries::RN)) {
            throw std::runtime_error(
                std::format("CramWriter: data series {} does not support name "
                            "tokeniser compression",
                            DataSeriesCode(ds)));
        }
        if ((method == CramBlockMethod::FQZCOMP) && !DataSeriesSupportsFqzcomp(ds)) {
            throw std::runtime_error(
                std::format("CramWriter: data series {} does not support fqzcomp compression",
                            DataSeriesCode(ds)));
        }

        return method;
    }

    CramBlockMethod NonDataSeriesMethod() const
    {
        switch (config.BlockCompressionMethod) {
            case CramBlockMethod::NAME_TOKENISER:
            case CramBlockMethod::FQZCOMP:
                return CramBlockMethod::GZIP;
            default:
                return config.BlockCompressionMethod;
        }
    }

    std::uint8_t EffectiveMinorVersion() const
    {
        std::uint8_t requiredMinor = 0;
        if (IsV31BlockMethod(NonDataSeriesMethod())) {
            requiredMinor = 1;
        }

        for (const auto& [series, _] : DATA_SERIES) {
            if (IsV31BlockMethod(ResolveDataSeriesMethod(series))) {
                requiredMinor = 1;
                break;
            }
        }
        return std::ranges::max(config.MinorVersion, requiredMinor);
    }

    void ValidateConfig() const
    {
        if (config.RecordsPerSlice <= 0) {
            throw std::runtime_error(std::format(
                "CramWriter: RecordsPerSlice must be positive, got {}", config.RecordsPerSlice));
        }
        if (config.CompressionLevel) {
            const int level = *config.CompressionLevel;
            if (level < 0 || level > 12) {
                throw std::runtime_error(
                    std::format("CramWriter: CompressionLevel must be in [0, 12], got {}", level));
            }
        }
    }

    void WaitForPendingWrite()
    {
        if (pendingWrite.valid()) {
            pendingWrite.get();
        }
    }

    void Open()
    {
        ValidateConfig();

        const int level = config.CompressionLevel.value_or(DEFAULT_GZIP_COMPRESSION_LEVEL);
        gzipCompressor.reset(libdeflate_alloc_compressor(level));
        if (!gzipCompressor) {
            throw std::runtime_error("CramWriter: failed to allocate gzip compressor");
        }

        writePath = detail::ResolveWritePath(path, config.UseTempFile);
        file.open(writePath, std::ios::binary | std::ios::trunc);
        if (!file) {
            throw std::runtime_error{std::format("CramWriter: cannot open {}", writePath.string())};
        }

        WriteFileDefinition();
        WriteHeaderContainer();
        nextContainerOffset = CurrentFileOffset();

        if (config.CompressionWorkers > 0) {
            compressionPool =
                std::make_shared<Parallel::ThreadPool<>>(Parallel::ThreadPool<>::Config{
                    .NumThreads = config.CompressionWorkers,
                });
        }

        pendingRecords.clear();
        pendingRecords.reserve(static_cast<std::size_t>(config.RecordsPerSlice));
        pendingSlices.clear();

        readGroupIndex.clear();
        for (std::int32_t i = 0; i < static_cast<std::int32_t>(std::size(header.ReadGroups()));
             ++i) {
            readGroupIndex.emplace(std::string{header.ReadGroups()[i].Id()}, i);
        }
    }

    void WriteRaw(std::span<const std::byte> data)
    {
        file.write(reinterpret_cast<const char*>(data.data()),
                   static_cast<std::streamsize>(std::size(data)));
        if (!file.good()) {
            throw std::runtime_error("CramWriter: write failed");
        }
    }

    void WriteSerializedContainer(std::vector<std::byte> data) { WriteRaw(data); }

    std::int64_t CurrentFileOffset()
    {
        const auto pos = file.tellp();
        if (pos == std::streampos{-1}) {
            throw std::runtime_error("CramWriter: failed to determine file offset");
        }
        return static_cast<std::int64_t>(pos);
    }

    std::filesystem::path CraiOutputPath() const
    {
        if (config.CraiPath) {
            return *config.CraiPath;
        }
        return DefaultCraiPath(path);
    }

    void AddCraiEntry(std::int32_t sequenceId, std::int64_t alignmentStart,
                      std::int64_t alignmentSpan, std::int64_t containerOffset,
                      std::int64_t sliceOffset, std::int64_t sliceSize)
    {
        craiEntries.push_back(CraiEntry{
            .SequenceId = sequenceId,
            .AlignmentStart = alignmentStart,
            .AlignmentSpan = alignmentSpan,
            .ContainerOffset = containerOffset,
            .SliceOffset = sliceOffset,
            .SliceSize = sliceSize,
        });
    }

    void CaptureCraiEntries(const std::vector<BamRecord>& sliceRecords,
                            const CramSliceHeader& sliceHeader, std::int64_t containerOffset,
                            std::int64_t sliceOffset, std::int64_t sliceSize)
    {
        if (!config.WriteCrai) {
            return;
        }

        const auto containerRefId = sliceHeader.RefSeqId;
        if (containerRefId >= 0) {
            const std::int64_t span = std::max<std::int64_t>(sliceHeader.AlignmentSpan, 1);
            AddCraiEntry(containerRefId, sliceHeader.AlignmentStart, span, containerOffset,
                         sliceOffset, sliceSize);
            return;
        }

        if (containerRefId == -1) {
            AddCraiEntry(-1, 0, 1, containerOffset, sliceOffset, sliceSize);
            return;
        }

        if (containerRefId != -2) {
            throw std::runtime_error(
                std::format("CramWriter: unsupported container ref ID {}", containerRefId));
        }

        std::vector<RefSpan> spans;
        std::unordered_map<std::int32_t, std::size_t> spanLookup;
        spans.reserve(std::size(sliceRecords));

        for (const BamRecord& record : sliceRecords) {
            if (!record.IsMapped() || record.RefId() < 0) {
                (void)LookupOrCreateRefSpan(spans, spanLookup, -1);
                continue;
            }

            RefSpan& span = LookupOrCreateRefSpan(spans, spanLookup, record.RefId());
            UpdateRefSpan(span, record);
        }

        for (const RefSpan& span : spans) {
            if (span.SequenceId == -1 || !span.HasMappedSpan) {
                AddCraiEntry(-1, 0, 1, containerOffset, sliceOffset, sliceSize);
                continue;
            }

            const std::int64_t alignmentStart = span.MinPos + 1;
            const std::int64_t alignmentSpan = std::max<std::int64_t>(span.MaxEnd - span.MinPos, 1);
            AddCraiEntry(span.SequenceId, alignmentStart, alignmentSpan, containerOffset,
                         sliceOffset, sliceSize);
        }
    }

    void WriteCraiSidecar()
    {
        if (!config.WriteCrai) {
            return;
        }

        std::string text;
        text.reserve(std::size(craiEntries) * 64);
        for (const CraiEntry& entry : craiEntries) {
            text += std::format("{}\t{}\t{}\t{}\t{}\t{}\n", entry.SequenceId, entry.AlignmentStart,
                                entry.AlignmentSpan, entry.ContainerOffset, entry.SliceOffset,
                                entry.SliceSize);
        }

        const std::span<const std::byte> payload{reinterpret_cast<const std::byte*>(text.data()),
                                                 std::size(text)};
        const std::vector<std::byte> compressed = CramGzipCompress(payload);

        const std::filesystem::path craiPath = CraiOutputPath();
        std::ofstream craiFile{craiPath, std::ios::binary | std::ios::trunc};
        if (!craiFile.is_open()) {
            throw std::runtime_error{
                std::format("CramWriter: cannot open CRAI sidecar {}", craiPath.string())};
        }

        craiFile.write(reinterpret_cast<const char*>(compressed.data()),
                       static_cast<std::streamsize>(std::size(compressed)));
        if (!craiFile.good()) {
            throw std::runtime_error{
                std::format("CramWriter: failed writing CRAI sidecar {}", craiPath.string())};
        }
    }

    void WriteFileDefinition()
    {
        CramFileDefinition def;
        def.MajorVersion = config.MajorVersion;
        def.MinorVersion = EffectiveMinorVersion();

        // Use a simple file ID (zeros are valid)
        const auto bytes = SerializeFileDefinition(def);
        WriteRaw(bytes);
    }

    void WriteHeaderContainer()
    {
        const auto headerText = header.ToText();

        // Build header block data: int32 length + text
        std::vector<std::byte> blockData;
        const auto textLen = static_cast<std::int32_t>(std::size(headerText));
        AppendLE(blockData, textLen);
        const auto* textPtr = reinterpret_cast<const std::byte*>(headerText.data());
        blockData.insert(std::end(blockData), textPtr, textPtr + textLen);

        // Build block
        CramBlock block;
        block.Method = CramBlockMethod::RAW;
        block.ContentType = CramBlockContentType::FILE_HEADER;
        block.ContentId = 0;
        block.RawSize = static_cast<std::int32_t>(std::size(blockData));
        block.CompressedSize = block.RawSize;
        block.Data = std::move(blockData);

        const auto blockBytes = SerializeBlock(block);

        // Build container header
        CramContainerHeader containerHeader;
        containerHeader.Length = static_cast<std::int32_t>(std::size(blockBytes));
        containerHeader.RefSeqId = 0;
        containerHeader.StartPos = 0;
        containerHeader.AlignmentSpan = 0;
        containerHeader.NumRecords = 0;
        containerHeader.RecordCounter = 0;
        containerHeader.Bases = 0;
        containerHeader.NumBlocks = 1;

        const auto containerBytes = SerializeContainerHeader(containerHeader);

        WriteRaw(containerBytes);
        WriteRaw(blockBytes);
    }

    struct PreparedSliceData
    {
        std::vector<RecordTagData> RecordTags;
    };

    struct ContainerEncodingPlan
    {
        CramCompressionHeader CompressionHeader;
        std::unordered_map<std::int32_t, CramBlockMethod> DataSeriesMethods;
        std::vector<PreparedSliceData> SliceData;
    };

    struct ContainerStats
    {
        std::int32_t RefSeqId{-2};
        std::int32_t StartPos{1};
        std::int32_t AlignmentSpan{0};
        std::int64_t Bases{0};
        std::int32_t NumRecords{0};
    };

    ContainerEncodingPlan BuildContainerEncodingPlan() const
    {
        ContainerEncodingPlan plan;
        plan.SliceData.resize(std::size(pendingSlices));

        plan.CompressionHeader.PreservationMap.ReadNamesIncluded = true;
        plan.CompressionHeader.PreservationMap.ApDelta = false;
        plan.CompressionHeader.PreservationMap.ReferenceRequired = false;

        plan.DataSeriesMethods.reserve(std::size(DATA_SERIES));
        for (const auto& [series, blockId] : DATA_SERIES) {
            plan.DataSeriesMethods.emplace(blockId, ResolveDataSeriesMethod(series));
        }

        plan.CompressionHeader.DataSeriesEncodings.reserve(std::size(DATA_SERIES));
        for (const auto& [series, blockId] : DATA_SERIES) {
            plan.CompressionHeader.DataSeriesEncodings.emplace_back(
                series, MakeDataSeriesEncodingDescriptor(series, blockId));
        }

        std::map<std::int32_t, char> tagContentTypes;
        std::vector<std::vector<TagTriple>> tagSets;
        tagSets.emplace_back();
        std::map<std::vector<TagTriple>, std::int32_t> tagSetIndex;
        tagSetIndex.emplace(std::vector<TagTriple>{}, 0);

        for (std::size_t sliceIndex = 0; sliceIndex < std::size(pendingSlices); ++sliceIndex) {
            const auto& sliceRecords = pendingSlices[sliceIndex];
            auto& prepared = plan.SliceData[sliceIndex];
            prepared.RecordTags.reserve(std::size(sliceRecords));

            for (const auto& record : sliceRecords) {
                auto tagData = BuildRecordTagData(record);
                for (const auto& payload : tagData.Payloads) {
                    if (!tagContentTypes.contains(payload.ContentId)) {
                        tagContentTypes.emplace(payload.ContentId, payload.Type);
                    }
                }

                if (const auto foundSet = tagSetIndex.find(tagData.TagSet);
                    foundSet != std::end(tagSetIndex)) {
                    tagData.TagListIndex = foundSet->second;
                } else {
                    const auto idx = static_cast<std::int32_t>(std::size(tagSets));
                    tagSetIndex.emplace(tagData.TagSet, idx);
                    tagSets.push_back(tagData.TagSet);
                    tagData.TagListIndex = idx;
                }

                if (const TagValue* rgValue = record.Tags().Get(RG_TAG)) {
                    if (const auto* rgText = std::get_if<std::string>(rgValue)) {
                        if (const auto rgIt = readGroupIndex.find(*rgText);
                            rgIt != std::end(readGroupIndex)) {
                            tagData.ReadGroupIndex = rgIt->second;
                        }
                    }
                }
                prepared.RecordTags.push_back(std::move(tagData));
            }
        }

        plan.CompressionHeader.PreservationMap.TagIdsDictionary = BuildTagIdsDictionary(tagSets);
        for (const auto& [contentId, type] : tagContentTypes) {
            plan.CompressionHeader.TagEncodings.emplace_back(
                contentId, MakeTagEncodingDescriptor(type, contentId));
        }

        return plan;
    }

    ContainerStats ComputeContainerStats() const
    {
        if (std::empty(pendingSlices)) {
            return {};
        }

        ContainerStats stats;
        std::int32_t minPos = std::numeric_limits<std::int32_t>::max();
        std::int32_t maxEnd = 0;
        bool haveFirstRef = false;
        std::int32_t firstRefId = -2;
        bool allSameRef = true;

        for (const auto& sliceRecords : pendingSlices) {
            for (const auto& rec : sliceRecords) {
                const auto refId = rec.RefId();
                if (!haveFirstRef) {
                    firstRefId = refId;
                    haveFirstRef = true;
                } else if (refId != firstRefId) {
                    allSameRef = false;
                }

                ++stats.NumRecords;
                stats.Bases += static_cast<std::int64_t>(std::size(rec.Sequence()));
                if (rec.Pos() >= 0) {
                    minPos = std::min(minPos, rec.Pos());
                    maxEnd = std::max(maxEnd, RecordAlignmentEnd(rec));
                }
            }
        }

        if (allSameRef && haveFirstRef) {
            stats.RefSeqId = firstRefId;
        } else {
            stats.RefSeqId = -2;
            minPos = 0;
            maxEnd = 0;
        }

        if (minPos == std::numeric_limits<std::int32_t>::max()) {
            minPos = 0;
        }
        stats.StartPos = minPos + 1;
        stats.AlignmentSpan = maxEnd - minPos;
        return stats;
    }

    CramSlice EncodeSlice(
        std::span<const BamRecord> records, std::span<const RecordTagData> recordTagData,
        const std::unordered_map<std::int32_t, CramBlockMethod>& dataSeriesMethods,
        std::int32_t containerRefId, std::int64_t recordCounterStart, bool useThreadLocalGzip,
        bool allowBlockParallel)
    {
        if (std::size(records) != std::size(recordTagData)) {
            throw std::runtime_error("CramWriter: per-slice record/tag metadata size mismatch");
        }

        const auto numRecords = static_cast<std::int32_t>(std::size(records));
        std::int32_t minPos = std::numeric_limits<std::int32_t>::max();
        std::int32_t maxEnd = 0;
        std::int64_t bases = 0;
        std::size_t totalNameBytes = 0;

        for (const auto& rec : records) {
            bases += static_cast<std::int64_t>(std::size(rec.Sequence()));
            totalNameBytes += std::size(rec.Name()) + 1;
            if (rec.Pos() >= 0) {
                minPos = std::min(minPos, rec.Pos());
                maxEnd = std::max(maxEnd, RecordAlignmentEnd(rec));
            }
        }
        if (containerRefId == -2) {
            minPos = 0;
            maxEnd = 0;
        } else if (minPos == std::numeric_limits<std::int32_t>::max()) {
            minPos = 0;
        }

        const bool needFqzMetadata =
            ResolveDataSeriesMethod(CramDataSeries::QS) == CramBlockMethod::FQZCOMP ||
            ResolveDataSeriesMethod(CramDataSeries::QQ) == CramBlockMethod::FQZCOMP;
        std::vector<std::uint32_t> qualityRecordLengths;
        std::vector<std::uint32_t> qualityRecordFlags;
        std::uint64_t totalQualityBytes = 0;
        if (needFqzMetadata) {
            qualityRecordLengths.reserve(std::size(records));
            qualityRecordFlags.reserve(std::size(records));
            for (const auto& rec : records) {
                const auto readLength = std::size(rec.Sequence());
                if (readLength > std::numeric_limits<std::uint32_t>::max()) {
                    throw std::runtime_error("CramWriter: read length exceeds fqz metadata limit");
                }
                qualityRecordLengths.push_back(static_cast<std::uint32_t>(readLength));
                totalQualityBytes += static_cast<std::uint32_t>(readLength);

                std::uint32_t flags = 0;
                if ((rec.Flag() & 0x10) != 0) {
                    flags |= FQZ_FREVERSE;
                }
                if ((rec.Flag() & 0x80) != 0) {
                    flags |= FQZ_FREAD2;
                }
                qualityRecordFlags.push_back(flags);
            }
        }

        CramExternalBlockStore extStore;
        const std::size_t recordReserve = static_cast<std::size_t>(numRecords) * 5U;
        extStore.ReserveBlock(BF_BLOCK_ID, recordReserve);
        extStore.ReserveBlock(CF_BLOCK_ID, recordReserve);
        extStore.ReserveBlock(RI_BLOCK_ID, recordReserve);
        extStore.ReserveBlock(RL_BLOCK_ID, recordReserve);
        extStore.ReserveBlock(AP_BLOCK_ID, recordReserve);
        extStore.ReserveBlock(RG_BLOCK_ID, recordReserve);
        extStore.ReserveBlock(RN_BLOCK_ID, totalNameBytes);
        extStore.ReserveBlock(MF_BLOCK_ID, recordReserve);
        extStore.ReserveBlock(NS_BLOCK_ID, recordReserve);
        extStore.ReserveBlock(NP_BLOCK_ID, recordReserve);
        extStore.ReserveBlock(TS_BLOCK_ID, recordReserve);
        extStore.ReserveBlock(TL_BLOCK_ID, recordReserve);
        extStore.ReserveBlock(FN_BLOCK_ID, recordReserve);
        extStore.ReserveBlock(FC_BLOCK_ID, static_cast<std::size_t>(numRecords) * 2U);
        extStore.ReserveBlock(FP_BLOCK_ID, recordReserve);
        extStore.ReserveBlock(MQ_BLOCK_ID, static_cast<std::size_t>(numRecords) * 2U);
        extStore.ReserveBlock(QS_BLOCK_ID, static_cast<std::size_t>(totalQualityBytes));
        if (bases > 0) {
            const std::size_t baseBytes = static_cast<std::size_t>(bases);
            extStore.ReserveBlock(BA_BLOCK_ID, baseBytes);
            extStore.ReserveBlock(BB_BLOCK_ID, baseBytes);
            extStore.ReserveBlock(IN_BLOCK_ID, baseBytes / 8U);
            extStore.ReserveBlock(SC_BLOCK_ID, baseBytes / 8U);
        }

        for (std::size_t i = 0; i < std::size(records); ++i) {
            EncodeRecord(records[i], recordTagData[i], extStore, containerRefId);
        }

        CramSliceHeader sliceHeader;
        sliceHeader.RefSeqId = containerRefId;
        sliceHeader.AlignmentStart = minPos + 1;
        sliceHeader.AlignmentSpan = maxEnd - minPos;
        sliceHeader.NumRecords = numRecords;
        sliceHeader.RecordCounter = recordCounterStart;
        sliceHeader.EmbeddedRefBlockId = -1;

        const auto contentIds = extStore.ContentIds();
        sliceHeader.BlockContentIds = contentIds;
        sliceHeader.NumBlocks = static_cast<std::int32_t>(std::size(contentIds)) + 1;

        std::vector<CramBlock> dataBlocks;
        std::vector<CramBlockMethod> dataBlockMethods;
        dataBlocks.reserve(1 + std::size(contentIds));
        dataBlockMethods.reserve(1 + std::size(contentIds));
        const auto nonDataSeriesMethod = NonDataSeriesMethod();

        CramBlock coreBlock;
        coreBlock.ContentType = CramBlockContentType::CORE_DATA;
        coreBlock.ContentId = 0;
        dataBlockMethods.push_back(ResolveBlockMethod(coreBlock.ContentType, coreBlock.ContentId,
                                                      dataSeriesMethods, nonDataSeriesMethod));
        dataBlocks.push_back(std::move(coreBlock));

        for (const auto id : contentIds) {
            CramBlock extBlock;
            extBlock.ContentType = CramBlockContentType::EXTERNAL_DATA;
            extBlock.ContentId = id;
            extBlock.Data = extStore.TakeBlockData(id);
            extBlock.RawSize = static_cast<std::int32_t>(std::size(extBlock.Data));
            dataBlockMethods.push_back(ResolveBlockMethod(extBlock.ContentType, extBlock.ContentId,
                                                          dataSeriesMethods, nonDataSeriesMethod));
            dataBlocks.push_back(std::move(extBlock));
        }

        if (std::size(dataBlocks) >
            static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max())) {
            throw std::runtime_error("CramWriter: too many blocks to compress in one slice");
        }

        const std::size_t workers = config.CompressionWorkers;
        constexpr std::size_t MIN_BLOCKS_FOR_PARALLEL = 8;
        const bool useParallelCompression = allowBlockParallel && compressionPool && workers > 1 &&
                                            std::size(dataBlocks) >= MIN_BLOCKS_FOR_PARALLEL;
        const bool useThreadLocalCompressor = useThreadLocalGzip || useParallelCompression;

        if (useParallelCompression) {
            std::atomic<std::int32_t> nextBlock{0};
            const std::int32_t totalBlocks = static_cast<std::int32_t>(std::size(dataBlocks));
            const std::int32_t taskCount =
                static_cast<std::int32_t>(std::ranges::min(workers, std::size(dataBlocks)));
            Parallel::Dispatch(
                compressionPool,
                MakeWorkStealingTask(
                    &nextBlock, totalBlocks,
                    [&dataBlocks, &dataBlockMethods, useThreadLocalCompressor,
                     &gzipCompressor = gzipCompressor, compressionLevel = config.CompressionLevel,
                     totalQualityBytes, qualityRecordLengths, qualityRecordFlags](std::int32_t i) {
                        CompressDataBlock(dataBlocks[static_cast<std::size_t>(i)],
                                          dataBlockMethods[static_cast<std::size_t>(i)],
                                          useThreadLocalCompressor, gzipCompressor,
                                          compressionLevel, totalQualityBytes, qualityRecordLengths,
                                          qualityRecordFlags);
                    }),
                taskCount);
        } else {
            for (std::int32_t i = 0; i < static_cast<std::int32_t>(std::size(dataBlocks)); ++i) {
                CompressDataBlock(dataBlocks[static_cast<std::size_t>(i)],
                                  dataBlockMethods[static_cast<std::size_t>(i)],
                                  useThreadLocalCompressor, gzipCompressor, config.CompressionLevel,
                                  totalQualityBytes, qualityRecordLengths, qualityRecordFlags);
            }
        }

        CramSlice slice;
        slice.Header = std::move(sliceHeader);
        bool sawCore = false;
        for (auto& block : dataBlocks) {
            if (block.ContentType == CramBlockContentType::CORE_DATA) {
                if (sawCore) {
                    throw std::runtime_error(
                        "CramWriter: duplicate CORE_DATA block in encoded slice");
                }
                slice.CoreBlock = std::move(block);
                sawCore = true;
            } else if (block.ContentType == CramBlockContentType::EXTERNAL_DATA) {
                slice.ExternalBlocks.push_back(std::move(block));
            } else {
                throw std::runtime_error("CramWriter: unexpected block content in encoded slice");
            }
        }
        if (!sawCore) {
            throw std::runtime_error("CramWriter: encoded slice missing CORE_DATA block");
        }

        if (slice.Header.NumBlocks !=
            static_cast<std::int32_t>(1 + std::size(slice.ExternalBlocks))) {
            throw std::runtime_error("CramWriter: encoded slice block count mismatch");
        }
        return slice;
    }

    void QueueContainerWrite(std::vector<std::byte> serializedContainer)
    {
        const auto bytesWritten = static_cast<std::int64_t>(std::size(serializedContainer));
        if (bytesWritten < 0 ||
            bytesWritten > (std::numeric_limits<std::int64_t>::max() - nextContainerOffset)) {
            throw std::runtime_error("CramWriter: container size overflow");
        }

        WaitForPendingWrite();
        if (compressionPool && config.CompressionWorkers > 0) {
            pendingWrite = std::async(std::launch::async, &Impl::WriteSerializedContainer, this,
                                      std::move(serializedContainer));
        } else {
            WriteRaw(serializedContainer);
        }
        nextContainerOffset += bytesWritten;
    }

    void FlushContainer()
    {
        if (std::empty(pendingSlices)) {
            return;
        }

        const ContainerStats containerStats = ComputeContainerStats();
        const ContainerEncodingPlan encodingPlan = BuildContainerEncodingPlan();

        std::vector<std::int64_t> sliceRecordCounters(std::size(pendingSlices), 0);
        std::int64_t runningRecordCounter = globalRecordCounter;
        for (std::size_t i = 0; i < std::size(pendingSlices); ++i) {
            sliceRecordCounters[i] = runningRecordCounter;
            runningRecordCounter += static_cast<std::int64_t>(std::size(pendingSlices[i]));
        }

        std::vector<CramSlice> slices(std::size(pendingSlices));
        const bool useParallelSlices =
            compressionPool && config.CompressionWorkers > 1 && std::size(pendingSlices) > 1;
        if (useParallelSlices) {
            std::atomic<std::size_t> nextSlice{0};
            const std::size_t sliceCount = std::size(pendingSlices);
            const auto taskCount = static_cast<std::int32_t>(
                std::ranges::min(config.CompressionWorkers, std::size(pendingSlices)));
            Parallel::Dispatch(
                compressionPool,
                MakeWorkStealingTask(
                    &nextSlice, sliceCount,
                    [this, &slices, &encodingPlan, &containerStats,
                     &sliceRecordCounters](std::size_t sliceIndex) {
                        slices[sliceIndex] =
                            EncodeSlice(pendingSlices[sliceIndex],
                                        encodingPlan.SliceData[sliceIndex].RecordTags,
                                        encodingPlan.DataSeriesMethods, containerStats.RefSeqId,
                                        sliceRecordCounters[sliceIndex],
                                        /*useThreadLocalGzip=*/true, /*allowBlockParallel=*/false);
                    }),
                taskCount);
        } else {
            for (std::size_t sliceIndex = 0; sliceIndex < std::size(pendingSlices); ++sliceIndex) {
                slices[sliceIndex] = EncodeSlice(
                    pendingSlices[sliceIndex], encodingPlan.SliceData[sliceIndex].RecordTags,
                    encodingPlan.DataSeriesMethods, containerStats.RefSeqId,
                    sliceRecordCounters[sliceIndex],
                    /*useThreadLocalGzip=*/false, /*allowBlockParallel=*/true);
            }
        }

        CramContainer container;
        container.Header.RefSeqId = containerStats.RefSeqId;
        container.Header.StartPos = containerStats.StartPos;
        container.Header.AlignmentSpan = containerStats.AlignmentSpan;
        container.Header.NumRecords = containerStats.NumRecords;
        container.Header.RecordCounter = globalRecordCounter;
        container.Header.Bases = containerStats.Bases;
        container.CompressionHeader = encodingPlan.CompressionHeader;
        container.Slices = std::move(slices);

        const std::int64_t containerOffset = nextContainerOffset;
        std::vector<std::byte> serialized = SerializeContainer(container);

        std::size_t containerHeaderBytesRead = 0;
        const auto serializedHeader = ParseContainerHeader(serialized, containerHeaderBytesRead);
        if (std::size(serializedHeader.Landmarks) != std::size(container.Slices)) {
            throw std::runtime_error("CramWriter: serialized landmark count mismatch");
        }

        for (std::size_t i = 0; i < std::size(container.Slices); ++i) {
            const auto sliceOffset = static_cast<std::int64_t>(serializedHeader.Landmarks[i]);
            std::int64_t sliceEnd = static_cast<std::int64_t>(serializedHeader.Length);
            if (i + 1 < std::size(serializedHeader.Landmarks)) {
                sliceEnd = static_cast<std::int64_t>(serializedHeader.Landmarks[i + 1]);
            }
            const auto sliceSize = sliceEnd - sliceOffset;
            if (sliceSize <= 0) {
                throw std::runtime_error("CramWriter: non-positive serialized slice size");
            }
            CaptureCraiEntries(pendingSlices[i], container.Slices[i].Header, containerOffset,
                               sliceOffset, sliceSize);
        }

        QueueContainerWrite(std::move(serialized));
        globalRecordCounter += containerStats.NumRecords;
        pendingSlices.clear();
    }

    void MaybeFlush()
    {
        if (static_cast<std::int32_t>(std::size(pendingRecords)) < config.RecordsPerSlice) {
            return;
        }

        pendingSlices.push_back(std::move(pendingRecords));
        pendingRecords.clear();
        pendingRecords.reserve(static_cast<std::size_t>(config.RecordsPerSlice));

        if (config.SlicesPerContainer > 0 &&
            static_cast<std::int32_t>(std::size(pendingSlices)) >= config.SlicesPerContainer) {
            FlushContainer();
        }
    }

    void EncodeRecord(const BamRecord& record, const RecordTagData& recordTagData,
                      CramExternalBlockStore& extStore, std::int32_t containerRefId)
    {
        const auto seq = record.Sequence();
        const auto readLength = static_cast<std::int32_t>(std::size(seq));
        const bool isUnmapped = (record.Flag() & 0x4) != 0;
        const auto qualities = record.Qualities();

        std::int32_t cramFlags = CRAM_FLAG_QUALITY_AS_ARRAY | CRAM_FLAG_DETACHED;
        if (std::empty(seq)) {
            cramFlags |= CRAM_FLAG_SEQUENCE_OMITTED;
        }

        extStore.WriteItf8(BF_BLOCK_ID, static_cast<std::int32_t>(record.Flag()));
        extStore.WriteItf8(CF_BLOCK_ID, cramFlags);
        if (containerRefId == -2) {
            extStore.WriteItf8(RI_BLOCK_ID, record.RefId());
        }
        extStore.WriteItf8(RL_BLOCK_ID, readLength);
        extStore.WriteItf8(AP_BLOCK_ID, record.Pos() + 1);
        extStore.WriteItf8(RG_BLOCK_ID, recordTagData.ReadGroupIndex);

        const auto name = record.Name();
        const auto* namePtr = reinterpret_cast<const std::byte*>(name.data());
        WriteByteArrayStop(extStore, RN_BLOCK_ID, {namePtr, std::size(name)}, std::byte{0});

        extStore.WriteItf8(MF_BLOCK_ID, 0);
        extStore.WriteItf8(NS_BLOCK_ID, record.NextRefId());
        extStore.WriteItf8(NP_BLOCK_ID, record.NextPos() + 1);
        extStore.WriteItf8(TS_BLOCK_ID, record.Tlen());
        extStore.WriteItf8(TL_BLOCK_ID, recordTagData.TagListIndex);

        for (const auto& payload : recordTagData.Payloads) {
            WriteEncodedTagPayload(extStore, payload);
        }

        if (!isUnmapped) {
            std::int32_t numFeatures = 0;
            std::int32_t seqOffset = 0;
            for (const auto& op : record.Cigar()) {
                const auto length = static_cast<std::int32_t>(op.Length());
                if (length <= 0) {
                    continue;
                }

                switch (op.Type()) {
                    case CigarOpType::M:
                    case CigarOpType::EQ:
                    case CigarOpType::X:
                    case CigarOpType::I:
                    case CigarOpType::S:
                        if (QueryCopyLength(readLength, length, seqOffset) > 0) {
                            ++numFeatures;
                        }
                        seqOffset += length;
                        break;
                    case CigarOpType::D:
                    case CigarOpType::N:
                    case CigarOpType::H:
                    case CigarOpType::P:
                        ++numFeatures;
                        break;
                }
            }
            if (numFeatures == 0 && !std::empty(seq)) {
                numFeatures = 1;
            }
            extStore.WriteItf8(FN_BLOCK_ID, numFeatures);

            std::int32_t readPos = 1;
            seqOffset = 0;
            std::int32_t previousFeaturePos = 0;
            bool wroteFeature = false;

            for (const auto& op : record.Cigar()) {
                const auto length = static_cast<std::int32_t>(op.Length());
                if (length <= 0) {
                    continue;
                }

                switch (op.Type()) {
                    case CigarOpType::M:
                    case CigarOpType::EQ:
                    case CigarOpType::X: {
                        const auto data = QuerySpan(seq, seqOffset, length);
                        if (!std::empty(data)) {
                            WriteFeatureHeader(extStore, 'b', readPos, previousFeaturePos);
                            WriteByteArrayLen(extStore, BB_BLOCK_ID, data);
                            wroteFeature = true;
                        }
                        readPos += length;
                        seqOffset += length;
                        break;
                    }
                    case CigarOpType::I: {
                        const auto data = QuerySpan(seq, seqOffset, length);
                        if (!std::empty(data)) {
                            WriteFeatureHeader(extStore, 'I', readPos, previousFeaturePos);
                            WriteByteArrayStop(extStore, IN_BLOCK_ID, data, std::byte{0});
                            wroteFeature = true;
                        }
                        readPos += length;
                        seqOffset += length;
                        break;
                    }
                    case CigarOpType::S: {
                        const auto data = QuerySpan(seq, seqOffset, length);
                        if (!std::empty(data)) {
                            WriteFeatureHeader(extStore, 'S', readPos, previousFeaturePos);
                            WriteByteArrayStop(extStore, SC_BLOCK_ID, data, std::byte{0});
                            wroteFeature = true;
                        }
                        readPos += length;
                        seqOffset += length;
                        break;
                    }
                    case CigarOpType::D:
                        WriteFeatureHeader(extStore, 'D', readPos, previousFeaturePos);
                        extStore.WriteItf8(DL_BLOCK_ID, length);
                        wroteFeature = true;
                        break;
                    case CigarOpType::N:
                        WriteFeatureHeader(extStore, 'N', readPos, previousFeaturePos);
                        extStore.WriteItf8(RS_BLOCK_ID, length);
                        wroteFeature = true;
                        break;
                    case CigarOpType::H:
                        WriteFeatureHeader(extStore, 'H', readPos, previousFeaturePos);
                        extStore.WriteItf8(HC_BLOCK_ID, length);
                        wroteFeature = true;
                        break;
                    case CigarOpType::P:
                        WriteFeatureHeader(extStore, 'P', readPos, previousFeaturePos);
                        extStore.WriteItf8(PD_BLOCK_ID, length);
                        wroteFeature = true;
                        break;
                }
            }

            if (!wroteFeature && !std::empty(seq)) {
                WriteFeatureHeader(extStore, 'b', 1, previousFeaturePos);
                const auto* seqPtr = reinterpret_cast<const std::byte*>(seq.data());
                WriteByteArrayLen(extStore, BB_BLOCK_ID, {seqPtr, std::size(seq)});
            }

            extStore.WriteItf8(MQ_BLOCK_ID, record.MapQ());
            WriteQualities(extStore, qualities, readLength);
            return;
        }

        const auto* seqPtr = reinterpret_cast<const std::byte*>(seq.data());
        extStore.WriteBytes(BA_BLOCK_ID, {seqPtr, std::size(seq)});
        WriteQualities(extStore, qualities, readLength);
    }

    void WriteEofContainer() { WriteRaw(CRAM_EOF_MARKER); }

    void Close()
    {
        if (closed) {
            return;
        }
        closed = true;

        if (!std::empty(pendingRecords)) {
            pendingSlices.push_back(std::move(pendingRecords));
            pendingRecords.clear();
        }

        FlushContainer();
        WaitForPendingWrite();

        WriteEofContainer();
        file.close();
        if (!file) {
            throw std::runtime_error(
                std::format("CramWriter: failed to close output {}", writePath.string()));
        }

        if (config.UseTempFile) {
            detail::AtomicRename(writePath, path, "CramWriter");
        }

        WriteCraiSidecar();
    }
};

// ---------------------------------------------------------------------------
// CramWriter
// ---------------------------------------------------------------------------

CramWriter::CramWriter(const std::filesystem::path& path, const SamHeader& header,
                       const CramWriterConfig& config)
    : impl_{std::make_unique<Impl>()}
{
    impl_->path = path;
    impl_->header = header;
    impl_->config = config;
    impl_->Open();
}

CramWriter::~CramWriter()
{
    if (impl_) {
        try {
            impl_->Close();
        } catch (...) {
            // suppress
        }
    }
}

CramWriter::CramWriter(CramWriter&&) noexcept = default;
CramWriter& CramWriter::operator=(CramWriter&&) noexcept = default;

void CramWriter::Write(const BamRecord& record)
{
    impl_->pendingRecords.push_back(record);
    impl_->MaybeFlush();
}

void CramWriter::Write(BamRecord&& record)
{
    impl_->pendingRecords.push_back(std::move(record));
    impl_->MaybeFlush();
}

void CramWriter::Write(const RawRecordView& record) { Write(record.ToOwned()); }

void CramWriter::Write(const RawRecord& record) { Write(record.View()); }

void CramWriter::WriteBatch(const RawRecordBatch& batch)
{
    for (std::size_t i{0}; i < batch.RecordCount(); ++i) {
        Write(batch.View(i));
    }
}

void CramWriter::Close() { impl_->Close(); }

}  // namespace Samoa
}  // namespace PacBio
