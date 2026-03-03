#include <pbsamoa/io/CramWriter.hpp>

#include <pbsamoa/cram/CramCodec.hpp>
#include <pbsamoa/cram/CramCompression.hpp>
#include <pbsamoa/cram/CramStructs.hpp>
#include <pbsamoa/index/CraiIndex.hpp>

#include "BinaryUtils.hpp"
#include "CramInternal.hpp"

#include <parallel/ThreadPool.h>

#include <htscodecs/fqzcomp_qual.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <format>
#include <fstream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <cstring>

namespace PacBio {
namespace Samoa {

namespace {

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

RecordTagData BuildRecordTagData(const BamRecord& record)
{
    RecordTagData result;

    for (const auto& [key, value] : record.Tags().Entries()) {
        if (key == RG_TAG) {
            continue;
        }

        const auto encoded = EncodeTagValueToBamPayload(value);
        const auto type = encoded.Type;
        result.TagSet.push_back(TagTriple{key.First(), key.Second(), type});

        EncodedTagPayload payload;
        payload.ContentId = TagContentId(key, type);
        payload.Type = type;
        payload.Data = encoded.Payload;
        result.Payloads.push_back(std::move(payload));
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
constexpr std::int32_t NF_BLOCK_ID = BlockIdFor(CramDataSeries::NF);
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

std::string DataSeriesCode(CramDataSeries ds)
{
    const auto value = std::to_underlying(ds);
    std::string code(2, '\0');
    code[0] = static_cast<char>((value >> 8) & 0xFF);
    code[1] = static_cast<char>(value & 0xFF);
    return code;
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

libdeflate_compressor* ThreadLocalGzipCompressor()
{
    const thread_local CramGzipCompressorContext THREAD_LOCAL_GZIP_COMPRESSOR_CONTEXT{};
    if (!THREAD_LOCAL_GZIP_COMPRESSOR_CONTEXT.Compressor) {
        throw std::runtime_error("CramWriter: failed to allocate gzip compressor");
    }
    return THREAD_LOCAL_GZIP_COMPRESSOR_CONTEXT.Compressor.get();
}

}  // namespace

// ---------------------------------------------------------------------------
// Impl
// ---------------------------------------------------------------------------

struct CramWriter::Impl
{
    std::filesystem::path path;
    std::ofstream file;
    SamHeader header;
    CramWriterConfig config;
    bool closed{false};

    std::vector<BamRecord> pendingRecords;
    std::int64_t globalRecordCounter{0};
    std::unordered_map<std::string, std::int32_t> readGroupIndex;
    std::vector<CraiEntry> craiEntries;
    CramGzipCompressorContext gzipContext;
    std::shared_ptr<Parallel::ThreadPool<>> compressionPool;

    CramBlockMethod DefaultDataSeriesMethod(CramDataSeries ds) const
    {
        switch (config.BlockCompressionMethod) {
            case CramBlockMethod::NAME_TOKENISER:
                return ds == CramDataSeries::RN ? CramBlockMethod::NAME_TOKENISER
                                                : CramBlockMethod::GZIP;
            case CramBlockMethod::FQZCOMP:
                return (ds == CramDataSeries::QS || ds == CramDataSeries::QQ)
                           ? CramBlockMethod::FQZCOMP
                           : CramBlockMethod::GZIP;
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

        if (method == CramBlockMethod::NAME_TOKENISER && ds != CramDataSeries::RN) {
            throw std::runtime_error(
                std::format("CramWriter: data series {} does not support name "
                            "tokeniser compression",
                            DataSeriesCode(ds)));
        }
        if (method == CramBlockMethod::FQZCOMP &&
            !(ds == CramDataSeries::QS || ds == CramDataSeries::QQ)) {
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

    void Open()
    {
        if (!gzipContext.Compressor) {
            throw std::runtime_error("CramWriter: failed to allocate gzip compressor");
        }

        file.open(path, std::ios::binary | std::ios::trunc);
        if (!file) {
            throw std::runtime_error{std::format("CramWriter: cannot open {}", path.string())};
        }

        WriteFileDefinition();
        WriteHeaderContainer();
        if (config.CompressionWorkers > 0) {
            compressionPool =
                std::make_shared<Parallel::ThreadPool<>>(Parallel::ThreadPool<>::Config{
                    .NumThreads = config.CompressionWorkers,
                });
        }

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
    }

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
        if (config.CraiPath.has_value()) {
            return *config.CraiPath;
        }
        return DefaultCraiPath(path);
    }

    static std::int64_t ReferenceEndForCrai(const BamRecord& record)
    {
        std::int64_t end = static_cast<std::int64_t>(record.Pos()) + 1;
        if (record.IsMapped()) {
            end = record.ReferenceEnd();
            if (end <= record.Pos()) {
                end = static_cast<std::int64_t>(record.Pos()) + 1;
            }
        }
        return end;
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

    void CaptureCraiEntries(std::int32_t containerRefId, const CramSliceHeader& sliceHeader,
                            std::int64_t containerOffset, std::int64_t sliceOffset,
                            std::int64_t sliceSize)
    {
        if (!config.WriteCrai) {
            return;
        }

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

        struct RefSpan
        {
            std::int32_t SequenceId{};
            bool HasMappedSpan{false};
            std::int64_t MinPos{};
            std::int64_t MaxEnd{};
        };

        std::vector<RefSpan> spans;
        std::unordered_map<std::int32_t, std::size_t> spanLookup;
        spans.reserve(std::size(pendingRecords));

        auto getSpan = [&](std::int32_t sequenceId) -> RefSpan& {
            if (const auto it = spanLookup.find(sequenceId); it != spanLookup.end()) {
                return spans[it->second];
            }
            const std::size_t idx = std::size(spans);
            spans.push_back(RefSpan{.SequenceId = sequenceId});
            spanLookup.emplace(sequenceId, idx);
            return spans.back();
        };

        for (const BamRecord& record : pendingRecords) {
            if (!record.IsMapped() || record.RefId() < 0) {
                (void)getSpan(-1);
                continue;
            }

            RefSpan& span = getSpan(record.RefId());
            const std::int64_t start = record.Pos();
            const std::int64_t end = ReferenceEndForCrai(record);
            if (!span.HasMappedSpan) {
                span.HasMappedSpan = true;
                span.MinPos = start;
                span.MaxEnd = end;
            } else {
                span.MinPos = std::min(span.MinPos, start);
                span.MaxEnd = std::max(span.MaxEnd, end);
            }
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
        WriteI32LE(blockData, textLen);
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

    void FlushPendingRecords()
    {
        if (std::empty(pendingRecords)) {
            return;
        }

        // Determine container properties
        std::int32_t containerRefId = -2;  // multi-ref by default
        std::int32_t minPos = std::numeric_limits<std::int32_t>::max();
        std::int32_t maxEnd = 0;
        std::int64_t bases = 0;
        std::size_t totalNameBytes = 0;

        bool allSameRef = true;
        const auto firstRefId = pendingRecords[0].RefId();
        for (const auto& rec : pendingRecords) {
            if (rec.RefId() != firstRefId) {
                allSameRef = false;
            }
            bases += static_cast<std::int64_t>(std::size(rec.Sequence()));
            totalNameBytes += std::size(rec.Name()) + 1;
            if (rec.Pos() >= 0) {
                if (rec.Pos() < minPos) {
                    minPos = rec.Pos();
                }
                std::int32_t end{rec.Pos() + static_cast<std::int32_t>(std::size(rec.Sequence()))};
                if (rec.IsMapped()) {
                    end = rec.ReferenceEnd();
                    if (end <= rec.Pos()) {
                        end = rec.Pos() + 1;
                    }
                }
                if (end > maxEnd) {
                    maxEnd = end;
                }
            }
        }

        if (allSameRef) {
            containerRefId = firstRefId;
        } else {
            minPos = 0;
            maxEnd = 0;
        }
        if (minPos == std::numeric_limits<std::int32_t>::max()) {
            minPos = 0;
        }

        const auto numRecords = static_cast<std::int32_t>(std::size(pendingRecords));
        const bool needFqzMetadata =
            ResolveDataSeriesMethod(CramDataSeries::QS) == CramBlockMethod::FQZCOMP ||
            ResolveDataSeriesMethod(CramDataSeries::QQ) == CramBlockMethod::FQZCOMP;

        std::vector<std::uint32_t> qualityRecordLengths;
        std::vector<std::uint32_t> qualityRecordFlags;
        std::uint64_t totalQualityBytes = 0;
        if (needFqzMetadata) {
            qualityRecordLengths.reserve(std::size(pendingRecords));
            qualityRecordFlags.reserve(std::size(pendingRecords));

            for (const auto& rec : pendingRecords) {
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

        // Build compression header with simple EXTERNAL encodings
        CramCompressionHeader compHeader;
        compHeader.PreservationMap.ReadNamesIncluded = true;
        compHeader.PreservationMap.ApDelta = false;
        compHeader.PreservationMap.ReferenceRequired = false;

        std::unordered_map<std::int32_t, CramBlockMethod> dataSeriesMethods;
        dataSeriesMethods.reserve(std::size(DATA_SERIES));
        for (const auto& [series, blockId] : DATA_SERIES) {
            dataSeriesMethods.emplace(blockId, ResolveDataSeriesMethod(series));
        }

        auto selectBlockMethod = [&](CramBlockContentType contentType,
                                     std::int32_t contentId) -> CramBlockMethod {
            if (contentType == CramBlockContentType::EXTERNAL_DATA) {
                if (const auto it = dataSeriesMethods.find(contentId);
                    it != std::end(dataSeriesMethods)) {
                    return it->second;
                }
            }
            return NonDataSeriesMethod();
        };

        auto makeExternalDesc = [](std::int32_t blockId) -> CramEncodingDescriptor {
            CramEncodingDescriptor desc;
            desc.CodecId = CramCodecId::EXTERNAL;
            WriteItf8(desc.Parameters, blockId);
            return desc;
        };
        auto makeByteArrayLenExternalDesc = [&](std::int32_t blockId) -> CramEncodingDescriptor {
            const auto lenDesc = makeExternalDesc(blockId);
            const auto dataDesc = makeExternalDesc(blockId);

            CramEncodingDescriptor desc;
            desc.CodecId = CramCodecId::BYTE_ARRAY_LEN;

            WriteItf8(desc.Parameters, std::to_underlying(lenDesc.CodecId));
            WriteItf8(desc.Parameters, static_cast<std::int32_t>(std::size(lenDesc.Parameters)));
            desc.Parameters.insert(std::end(desc.Parameters), std::begin(lenDesc.Parameters),
                                   std::end(lenDesc.Parameters));

            WriteItf8(desc.Parameters, std::to_underlying(dataDesc.CodecId));
            WriteItf8(desc.Parameters, static_cast<std::int32_t>(std::size(dataDesc.Parameters)));
            desc.Parameters.insert(std::end(desc.Parameters), std::begin(dataDesc.Parameters),
                                   std::end(dataDesc.Parameters));

            return desc;
        };
        auto makeByteArrayStopExternalDesc = [&](std::int32_t blockId,
                                                 std::byte stopByte) -> CramEncodingDescriptor {
            CramEncodingDescriptor desc;
            desc.CodecId = CramCodecId::BYTE_ARRAY_STOP;
            desc.Parameters.push_back(stopByte);
            WriteItf8(desc.Parameters, blockId);
            return desc;
        };
        auto makeTagDesc = [&](char type, std::int32_t blockId) -> CramEncodingDescriptor {
            if (type == 'Z' || type == 'H') {
                return makeByteArrayStopExternalDesc(blockId, std::byte{'\t'});
            }
            return makeByteArrayLenExternalDesc(blockId);
        };

        auto makeDataSeriesEncoding = [&](CramDataSeries series,
                                          std::int32_t blockId) -> CramEncodingDescriptor {
            switch (series) {
                case CramDataSeries::RN:
                case CramDataSeries::IN:
                case CramDataSeries::SC:
                    return makeByteArrayStopExternalDesc(blockId, std::byte{0});
                case CramDataSeries::BB:
                    return makeByteArrayLenExternalDesc(blockId);
                default:
                    return makeExternalDesc(blockId);
            }
        };

        compHeader.DataSeriesEncodings.reserve(std::size(DATA_SERIES));
        for (const auto& [series, blockId] : DATA_SERIES) {
            compHeader.DataSeriesEncodings.emplace_back(series,
                                                        makeDataSeriesEncoding(series, blockId));
        }

        // Build per-record tag payloads and tag-set dictionary.
        std::vector<RecordTagData> recordTagData;
        recordTagData.reserve(std::size(pendingRecords));

        std::map<std::int32_t, char> tagContentTypes;
        std::vector<std::vector<TagTriple>> tagSets;
        tagSets.emplace_back();  // index 0: empty tag set

        std::map<std::vector<TagTriple>, std::int32_t> tagSetIndex;
        tagSetIndex.emplace(std::vector<TagTriple>{}, 0);

        for (const auto& record : pendingRecords) {
            auto tagData = BuildRecordTagData(record);
            for (const auto& payload : tagData.Payloads) {
                if (!tagContentTypes.contains(payload.ContentId)) {
                    tagContentTypes.emplace(payload.ContentId, payload.Type);
                }
            }
            const auto foundSet = tagSetIndex.find(tagData.TagSet);
            if (foundSet == std::end(tagSetIndex)) {
                const auto idx = static_cast<std::int32_t>(std::size(tagSets));
                tagSetIndex.emplace(tagData.TagSet, idx);
                tagSets.push_back(tagData.TagSet);
                tagData.TagListIndex = idx;
            } else {
                tagData.TagListIndex = foundSet->second;
            }

            if (const TagValue* rgValue = record.Tags().Get(RG_TAG)) {
                if (const auto* rgText = std::get_if<std::string>(rgValue)) {
                    if (const auto rgIt = readGroupIndex.find(*rgText);
                        rgIt != std::end(readGroupIndex)) {
                        tagData.ReadGroupIndex = rgIt->second;
                    }
                }
            }
            recordTagData.push_back(std::move(tagData));
        }

        compHeader.PreservationMap.TagIdsDictionary = BuildTagIdsDictionary(tagSets);

        for (const auto& [contentId, type] : tagContentTypes) {
            compHeader.TagEncodings.emplace_back(contentId, makeTagDesc(type, contentId));
        }

        // Encode records into blocks
        CramBitWriter coreWriter;  // unused with all-external encoding
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

        for (std::size_t i = 0; i < std::size(pendingRecords); ++i) {
            EncodeRecord(pendingRecords[i], recordTagData[i], extStore, containerRefId);
        }
        coreWriter.Flush();

        // Serialize compression header block
        const auto compHeaderData = SerializeCompressionHeader(compHeader);
        CramBlock compHdrBlock;
        compHdrBlock.Method = CramBlockMethod::RAW;
        compHdrBlock.ContentType = CramBlockContentType::COMPRESSION_HEADER;
        compHdrBlock.ContentId = 0;
        compHdrBlock.RawSize = static_cast<std::int32_t>(std::size(compHeaderData));
        compHdrBlock.CompressedSize = compHdrBlock.RawSize;
        compHdrBlock.Data = compHeaderData;
        const auto compHdrBytes = SerializeBlock(compHdrBlock);

        // Build slice header
        CramSliceHeader sliceHeader;
        sliceHeader.RefSeqId = containerRefId;
        sliceHeader.AlignmentStart = minPos + 1;  // 1-based
        sliceHeader.AlignmentSpan = maxEnd - minPos;
        sliceHeader.NumRecords = numRecords;
        sliceHeader.RecordCounter = globalRecordCounter;
        sliceHeader.EmbeddedRefBlockId = -1;

        // Collect external block IDs
        const auto contentIds = extStore.ContentIds();
        sliceHeader.BlockContentIds = contentIds;
        // +1 for core data block
        sliceHeader.NumBlocks = static_cast<std::int32_t>(std::size(contentIds)) + 1;

        // Serialize slice header block
        const auto sliceHdrData = SerializeSliceHeader(sliceHeader);
        CramBlock sliceHdrBlock;
        sliceHdrBlock.Method = CramBlockMethod::RAW;
        sliceHdrBlock.ContentType = CramBlockContentType::SLICE_HEADER;
        sliceHdrBlock.ContentId = 0;
        sliceHdrBlock.RawSize = static_cast<std::int32_t>(std::size(sliceHdrData));
        sliceHdrBlock.CompressedSize = sliceHdrBlock.RawSize;
        sliceHdrBlock.Data = sliceHdrData;
        const auto sliceHdrBytes = SerializeBlock(sliceHdrBlock);

        std::vector<CramBlock> dataBlocks;
        std::vector<CramBlockMethod> dataBlockMethods;
        dataBlocks.reserve(1 + std::size(contentIds));
        dataBlockMethods.reserve(1 + std::size(contentIds));

        // Core data block
        CramBlock coreBlock;
        coreBlock.ContentType = CramBlockContentType::CORE_DATA;
        coreBlock.ContentId = 0;
        coreBlock.Data = std::move(coreWriter).Data();
        coreBlock.RawSize = static_cast<std::int32_t>(std::size(coreBlock.Data));
        dataBlockMethods.push_back(selectBlockMethod(coreBlock.ContentType, coreBlock.ContentId));
        dataBlocks.push_back(std::move(coreBlock));

        // External data blocks
        for (const auto id : contentIds) {
            CramBlock extBlock;
            extBlock.ContentType = CramBlockContentType::EXTERNAL_DATA;
            extBlock.ContentId = id;
            extBlock.Data = extStore.TakeBlockData(id);
            extBlock.RawSize = static_cast<std::int32_t>(std::size(extBlock.Data));
            dataBlockMethods.push_back(selectBlockMethod(extBlock.ContentType, extBlock.ContentId));
            dataBlocks.push_back(std::move(extBlock));
        }

        if (std::size(dataBlocks) >
            static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max())) {
            throw std::runtime_error("CramWriter: too many blocks to compress in one slice");
        }

        const std::size_t workers = config.CompressionWorkers;
        constexpr std::size_t MIN_BLOCKS_FOR_PARALLEL = 8;
        const bool useParallelCompression = static_cast<bool>(compressionPool) && workers > 1 &&
                                            std::size(dataBlocks) >= MIN_BLOCKS_FOR_PARALLEL;
        auto compressDataBlock = [&](std::int32_t blockIndex) {
            auto& block = dataBlocks[static_cast<std::size_t>(blockIndex)];
            const auto method = dataBlockMethods[static_cast<std::size_t>(blockIndex)];

            if (method == CramBlockMethod::FQZCOMP &&
                (block.ContentId == QS_BLOCK_ID || block.ContentId == QQ_BLOCK_ID)) {
                if (totalQualityBytes != static_cast<std::uint64_t>(block.RawSize)) {
                    throw std::runtime_error(
                        std::format("CramWriter: fqz metadata size mismatch for content ID {} "
                                    "(metadata={}, raw={})",
                                    block.ContentId, totalQualityBytes, block.RawSize));
                }
                block.Data =
                    CramFqzcompCompress(block.Data, qualityRecordLengths, qualityRecordFlags);
                block.Method = CramBlockMethod::FQZCOMP;
                if (std::size(block.Data) > std::numeric_limits<std::int32_t>::max()) {
                    throw std::runtime_error(
                        "CramWriter: fqz compressed block exceeds CRAM size "
                        "limit");
                }
                block.CompressedSize = static_cast<std::int32_t>(std::size(block.Data));
                return;
            }

            CompressCramBlock(block, method,
                              useParallelCompression ? ThreadLocalGzipCompressor()
                                                     : gzipContext.Compressor.get());
        };

        if (useParallelCompression) {
            std::atomic<std::int32_t> nextBlock{0};
            const std::int32_t totalBlocks = static_cast<std::int32_t>(std::size(dataBlocks));
            const std::int32_t taskCount =
                static_cast<std::int32_t>(std::ranges::min(workers, std::size(dataBlocks)));
            Parallel::Dispatch(
                compressionPool,
                [&, totalBlocks](std::int32_t) {
                    while (true) {
                        const std::int32_t blockIndex =
                            nextBlock.fetch_add(1, std::memory_order_relaxed);
                        if (blockIndex >= totalBlocks) {
                            break;
                        }
                        compressDataBlock(blockIndex);
                    }
                },
                taskCount);
        } else {
            for (std::int32_t i = 0; i < static_cast<std::int32_t>(std::size(dataBlocks)); ++i) {
                compressDataBlock(i);
            }
        }

        std::vector<std::vector<std::byte>> dataBlockBytes;
        dataBlockBytes.reserve(std::size(dataBlocks));
        for (auto& block : dataBlocks) {
            dataBlockBytes.push_back(SerializeBlock(block));
        }

        // Calculate total blocks data size
        std::int32_t totalBlocksSize = static_cast<std::int32_t>(std::size(compHdrBytes));
        totalBlocksSize += static_cast<std::int32_t>(std::size(sliceHdrBytes));
        for (const auto& blockBytes : dataBlockBytes) {
            totalBlocksSize += static_cast<std::int32_t>(std::size(blockBytes));
        }

        // Container header
        CramContainerHeader containerHdr;
        containerHdr.Length = totalBlocksSize;
        containerHdr.RefSeqId = containerRefId;
        containerHdr.StartPos = minPos + 1;  // 1-based
        containerHdr.AlignmentSpan = maxEnd - minPos;
        containerHdr.NumRecords = numRecords;
        containerHdr.RecordCounter = globalRecordCounter;
        containerHdr.Bases = bases;
        containerHdr.NumBlocks = 3 + static_cast<std::int32_t>(std::size(
                                         contentIds));  // comp hdr + slice hdr + core + externals
        // Landmark: offset of slice from end of container header
        containerHdr.Landmarks = {static_cast<std::int32_t>(std::size(compHdrBytes))};

        const auto containerHdrBytes = SerializeContainerHeader(containerHdr);
        const std::int64_t containerOffset = CurrentFileOffset();
        const std::int64_t sliceOffset = containerHdr.Landmarks.at(0);
        const std::int64_t sliceSize = static_cast<std::int64_t>(totalBlocksSize) -
                                       static_cast<std::int64_t>(std::size(compHdrBytes));

        // Write everything
        WriteRaw(containerHdrBytes);
        WriteRaw(compHdrBytes);
        WriteRaw(sliceHdrBytes);
        for (const auto& blockBytes : dataBlockBytes) {
            WriteRaw(blockBytes);
        }

        CaptureCraiEntries(containerRefId, sliceHeader, containerOffset, sliceOffset, sliceSize);

        globalRecordCounter += numRecords;
        pendingRecords.clear();
    }

    void MaybeFlush()
    {
        if (static_cast<std::int32_t>(std::size(pendingRecords)) >= config.RecordsPerSlice) {
            FlushPendingRecords();
        }
    }

    void EncodeRecord(const BamRecord& record, const RecordTagData& recordTagData,
                      CramExternalBlockStore& extStore, std::int32_t containerRefId)
    {
        const auto seq = record.Sequence();
        const auto readLength = static_cast<std::int32_t>(std::size(seq));
        const bool isUnmapped = (record.Flag() & 0x4) != 0;

        auto queryCopyLength = [readLength](std::int32_t requested,
                                            std::int32_t seqOffset) -> std::int32_t {
            if (requested <= 0 || seqOffset < 0 || seqOffset >= readLength) {
                return 0;
            }
            return std::ranges::min(requested, readLength - seqOffset);
        };
        auto querySpan = [&](std::int32_t seqOffset,
                             std::int32_t requested) -> std::span<const std::byte> {
            const auto copyLen = queryCopyLength(requested, seqOffset);
            if (copyLen <= 0) {
                return {};
            }
            const auto* begin = reinterpret_cast<const std::byte*>(seq.data() + seqOffset);
            return {begin, static_cast<std::size_t>(copyLen)};
        };
        auto writeByteArrayLen = [&](std::int32_t blockId, std::span<const std::byte> data) {
            extStore.WriteItf8(blockId, static_cast<std::int32_t>(std::size(data)));
            extStore.WriteBytes(blockId, data);
        };
        auto writeByteArrayStop = [&](std::int32_t blockId, std::span<const std::byte> data,
                                      std::byte stop) {
            extStore.WriteBytes(blockId, data);
            extStore.WriteByte(blockId, stop);
        };
        auto writeQualities = [&](std::int32_t expectedLength) {
            const auto quals = record.Qualities();
            if (static_cast<std::int32_t>(std::size(quals)) == expectedLength) {
                const auto* qPtr = reinterpret_cast<const std::byte*>(quals.data());
                extStore.WriteBytes(QS_BLOCK_ID, {qPtr, std::size(quals)});
                return;
            }
            for (std::int32_t qi = 0; qi < expectedLength; ++qi) {
                const auto q = (qi < static_cast<std::int32_t>(std::size(quals)))
                                   ? quals[qi]
                                   : static_cast<std::uint8_t>(0xFF);
                extStore.WriteByte(QS_BLOCK_ID, static_cast<std::byte>(q));
            }
        };

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
        writeByteArrayStop(RN_BLOCK_ID, {namePtr, std::size(name)}, std::byte{0});

        extStore.WriteItf8(MF_BLOCK_ID, 0);
        extStore.WriteItf8(NS_BLOCK_ID, record.NextRefId());
        extStore.WriteItf8(NP_BLOCK_ID, record.NextPos() + 1);
        extStore.WriteItf8(TS_BLOCK_ID, record.Tlen());
        extStore.WriteItf8(TL_BLOCK_ID, recordTagData.TagListIndex);

        for (const auto& payload : recordTagData.Payloads) {
            if (payload.Type == 'Z' || payload.Type == 'H') {
                writeByteArrayStop(payload.ContentId, payload.Data, std::byte{'\t'});
            } else {
                writeByteArrayLen(payload.ContentId, payload.Data);
            }
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
                        if (queryCopyLength(length, seqOffset) > 0) {
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

            auto writeFeatureHeader = [&](char code, std::int32_t featurePos) {
                extStore.WriteByte(FC_BLOCK_ID,
                                   static_cast<std::byte>(static_cast<std::uint8_t>(code)));
                extStore.WriteItf8(FP_BLOCK_ID, featurePos - previousFeaturePos);
                previousFeaturePos = featurePos;
            };

            for (const auto& op : record.Cigar()) {
                const auto length = static_cast<std::int32_t>(op.Length());
                if (length <= 0) {
                    continue;
                }

                switch (op.Type()) {
                    case CigarOpType::M:
                    case CigarOpType::EQ:
                    case CigarOpType::X: {
                        const auto data = querySpan(seqOffset, length);
                        if (!std::empty(data)) {
                            writeFeatureHeader('b', readPos);
                            writeByteArrayLen(BB_BLOCK_ID, data);
                            wroteFeature = true;
                        }
                        readPos += length;
                        seqOffset += length;
                        break;
                    }
                    case CigarOpType::I: {
                        const auto data = querySpan(seqOffset, length);
                        if (!std::empty(data)) {
                            writeFeatureHeader('I', readPos);
                            writeByteArrayStop(IN_BLOCK_ID, data, std::byte{0});
                            wroteFeature = true;
                        }
                        readPos += length;
                        seqOffset += length;
                        break;
                    }
                    case CigarOpType::S: {
                        const auto data = querySpan(seqOffset, length);
                        if (!std::empty(data)) {
                            writeFeatureHeader('S', readPos);
                            writeByteArrayStop(SC_BLOCK_ID, data, std::byte{0});
                            wroteFeature = true;
                        }
                        readPos += length;
                        seqOffset += length;
                        break;
                    }
                    case CigarOpType::D:
                        writeFeatureHeader('D', readPos);
                        extStore.WriteItf8(DL_BLOCK_ID, length);
                        wroteFeature = true;
                        break;
                    case CigarOpType::N:
                        writeFeatureHeader('N', readPos);
                        extStore.WriteItf8(RS_BLOCK_ID, length);
                        wroteFeature = true;
                        break;
                    case CigarOpType::H:
                        writeFeatureHeader('H', readPos);
                        extStore.WriteItf8(HC_BLOCK_ID, length);
                        wroteFeature = true;
                        break;
                    case CigarOpType::P:
                        writeFeatureHeader('P', readPos);
                        extStore.WriteItf8(PD_BLOCK_ID, length);
                        wroteFeature = true;
                        break;
                }
            }

            if (!wroteFeature && !std::empty(seq)) {
                writeFeatureHeader('b', 1);
                const auto* seqPtr = reinterpret_cast<const std::byte*>(seq.data());
                writeByteArrayLen(BB_BLOCK_ID, {seqPtr, std::size(seq)});
            }

            extStore.WriteItf8(MQ_BLOCK_ID, record.MapQ());
            writeQualities(readLength);
            return;
        }

        const auto* seqPtr = reinterpret_cast<const std::byte*>(seq.data());
        extStore.WriteBytes(BA_BLOCK_ID, {seqPtr, std::size(seq)});
        writeQualities(readLength);
    }

    void WriteEofContainer() { WriteRaw(CRAM_EOF_MARKER); }

    void Close()
    {
        if (closed) {
            return;
        }
        closed = true;

        FlushPendingRecords();

        WriteEofContainer();
        file.close();
        if (!file) {
            throw std::runtime_error(
                std::format("CramWriter: failed to close output {}", path.string()));
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

void CramWriter::Write(const RawRecord& record) { Write(record.ToOwned()); }

void CramWriter::WriteBatch(const RawRecordBatch& batch)
{
    for (std::size_t i{0}; i < batch.RecordCount(); ++i) {
        Write(RawRecord{batch.RecordData(i)});
    }
}

void CramWriter::Close() { impl_->Close(); }

}  // namespace Samoa
}  // namespace PacBio
