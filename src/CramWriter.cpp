#include <pbsamoa/io/CramWriter.hpp>

#include <pbsamoa/cram/CramCodec.hpp>
#include <pbsamoa/cram/CramCompression.hpp>
#include <pbsamoa/cram/CramStructs.hpp>
#include <pbsamoa/index/CraiIndex.hpp>

#include "BinaryUtils.hpp"
#include "CramInternal.hpp"
#include "CramStructs.hpp"
#include "Md5.hpp"
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
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

#include <cstring>

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
constexpr std::int32_t BS_BLOCK_ID = BlockIdFor(CramDataSeries::BS);
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
    if ((method == CramBlockMethod::RANS4X16) &&
        ((block.ContentId == QS_BLOCK_ID) || (block.ContentId == QQ_BLOCK_ID))) {
        CompressRansQualityBlock(block);
        return;
    }
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

// The code of a substitution of the reference base by the read base in the default substitution
// matrix (each row lists the other bases of ACGTN in order: code 0..3), or -1 when one of the
// bases is not in ACGTN.
std::int32_t SubstitutionCode(const char refBase, const char readBase)
{
    constexpr std::string_view BASES{"ACGTN"};
    const std::size_t ref{BASES.find(refBase)};
    const std::size_t read{BASES.find(readBase)};
    if ((ref == std::string_view::npos) || (read == std::string_view::npos) || (ref == read)) {
        return -1;
    }
    return static_cast<std::int32_t>((read < ref) ? read : (read - 1));
}

char UppercaseBase(const char base)
{
    return ((base >= 'a') && (base <= 'z')) ? static_cast<char>(base - ('a' - 'A')) : base;
}

// The reference of CramWriterConfig::ReferenceSequences: views of the text, no copy.
class TextReference final : public CramReference
{
public:
    explicit TextReference(std::vector<std::string_view> sequences)
        : sequences_{std::move(sequences)}
    {
    }

    std::int32_t NumSequences() const override
    {
        return static_cast<std::int32_t>(std::size(sequences_));
    }

    std::int64_t Length(const std::int32_t refId) const override
    {
        return static_cast<std::int64_t>(std::size(sequences_[static_cast<std::size_t>(refId)]));
    }

    std::string_view Bases(const std::int32_t refId, const std::int64_t start,
                           const std::int64_t length, std::string& /*buffer*/) const override
    {
        return sequences_[static_cast<std::size_t>(refId)].substr(static_cast<std::size_t>(start),
                                                                  static_cast<std::size_t>(length));
    }

private:
    std::vector<std::string_view> sequences_;
};

// The MD5 of the bases [start, start + length) of a reference in upper case, from windows of
// 64 KiB: a slice can span megabases, and a packed reference decodes each window to buffer.
std::array<std::byte, 16> ReferenceMd5(const CramReference& reference, const std::int32_t refId,
                                       const std::int64_t start, const std::int64_t length,
                                       std::string& buffer)
{
    constexpr std::int64_t WINDOW_BASES{std::int64_t{1} << 16};
    const std::int64_t end{start + length};
    std::int64_t next{start};
    return detail::ComputeUppercaseMd5([&]() -> std::string_view {
        if (next >= end) {
            return {};
        }
        const std::int64_t n{std::min(WINDOW_BASES, end - next)};
        const std::string_view bases{reference.Bases(refId, next, n, buffer)};
        next += n;
        return bases;
    });
}

// The reference bases that the CIGAR operations of a record cover.
std::int64_t ReferenceSpan(const BamRecord& record)
{
    std::int64_t span{0};
    for (const auto& op : record.Cigar()) {
        switch (op.Type()) {
            case CigarOpType::M:
            case CigarOpType::EQ:
            case CigarOpType::X:
            case CigarOpType::D:
            case CigarOpType::N:
                span += op.Length();
                break;
            default:
                break;
        }
    }
    return span;
}

namespace {

constexpr std::int32_t REFERENCE_MATCH_BYTES{8};

bool ReferenceWordMatches(const std::span<const std::byte> bases, const std::string_view reference,
                          const std::int64_t refPos)
{
    if ((std::ssize(bases) < REFERENCE_MATCH_BYTES) || (refPos < 0) ||
        (refPos > (std::ssize(reference) - REFERENCE_MATCH_BYTES))) {
        return false;
    }
    std::uint64_t readWord{0};
    std::uint64_t referenceWord{0};
    std::memcpy(&readWord, bases.data(), REFERENCE_MATCH_BYTES);
    std::memcpy(&referenceWord, reference.data() + refPos, REFERENCE_MATCH_BYTES);
    // Lowercase reference bases require the scalar case conversion below.
    constexpr std::uint64_t LOWERCASE_MASK{0x2020202020202020ULL};
    return ((referenceWord & LOWERCASE_MASK) == 0) && (readWord == referenceWord);
}

}  // namespace

// Writes the read features of the bases of an M, = or X operation against the reference: nothing
// for a base equal to the reference, a substitution (X) where both bases are in ACGTN, else a run
// of bases (b), e.g. for a base past the end of the reference or an IUPAC code.
void WriteBasesAgainstReference(CramExternalBlockStore& extStore, std::span<const std::byte> bases,
                                const std::string_view reference, const std::int64_t refPos,
                                const std::int32_t readPos, std::int32_t& previousFeaturePos,
                                std::int32_t& numFeatures)
{
    std::int32_t runStart{-1};
    const auto flushRun = [&](const std::int32_t end) {
        if (runStart >= 0) {
            WriteFeatureHeader(extStore, 'b', readPos + runStart, previousFeaturePos);
            WriteByteArrayLen(extStore, BB_BLOCK_ID,
                              bases.subspan(static_cast<std::size_t>(runStart),
                                            static_cast<std::size_t>(end - runStart)));
            ++numFeatures;
            runStart = -1;
        }
    };
    const auto n = static_cast<std::int32_t>(std::size(bases));
    for (std::int32_t i = 0; i < n; ++i) {
        const std::int64_t r{refPos + i};
        if (ReferenceWordMatches(bases.subspan(i), reference, r)) {
            flushRun(i);
            i += REFERENCE_MATCH_BYTES - 1;
            continue;
        }
        const auto readBase = static_cast<char>(std::to_integer<std::uint8_t>(bases[i]));
        const bool onReference{(r >= 0) && (r < static_cast<std::int64_t>(std::size(reference)))};
        const char refBase{onReference ? UppercaseBase(reference[static_cast<std::size_t>(r)])
                                       : '\0'};
        if (onReference && (readBase == refBase)) {
            flushRun(i);
            continue;
        }
        const std::int32_t code{onReference ? SubstitutionCode(refBase, readBase) : -1};
        if (code >= 0) {
            flushRun(i);
            WriteFeatureHeader(extStore, 'X', readPos + i, previousFeaturePos);
            extStore.WriteByte(BS_BLOCK_ID, static_cast<std::byte>(code));
            ++numFeatures;
            continue;
        }
        if (runStart < 0) {
            runStart = i;
        }
    }
    flushRun(n);
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
    // the reference of reference-based CRAM (Reference, or a view of ReferenceSequences)
    std::shared_ptr<const CramReference> reference;
    bool closed{false};

    std::vector<BamRecord> pendingRecords;
    std::vector<std::vector<BamRecord>> pendingSlices;
    std::int64_t globalRecordCounter{0};
    std::int64_t nextContainerOffset{0};
    std::unordered_map<std::string, std::int32_t> readGroupIndex;
    std::vector<CraiEntry> craiEntries;
    std::future<void> pendingWrite;
    // the compressor of the writer thread; EncodeContainer uses thread-local ones
    mutable LibdeflateCompressorPtr gzipCompressor;
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

    bool ReferenceBased() const { return reference != nullptr; }

    // the reference of a mapped record with a reference sequence of at least one base
    bool HasReference(const BamRecord& record) const
    {
        const std::int32_t refId{record.RefId()};
        return ReferenceBased() && (refId >= 0) && (refId < reference->NumSequences()) &&
               (reference->Length(refId) > 0);
    }

    // Reference-based CRAM: one reference per @SQ line, of its length; M5 where it is missing.
    void PrepareReferences()
    {
        if (!std::empty(config.ReferenceSequences) && config.Reference) {
            throw std::runtime_error{"CramWriter: set ReferenceSequences or Reference, not both"};
        }
        reference = config.Reference;
        if (!std::empty(config.ReferenceSequences)) {
            reference = std::make_shared<const TextReference>(config.ReferenceSequences);
        }
        if (!ReferenceBased()) {
            return;
        }
        std::vector<ReferenceSequence>& sequences{header.ReferenceSequences()};
        if (std::cmp_not_equal(std::size(sequences), reference->NumSequences())) {
            throw std::runtime_error{
                std::format("CramWriter: {} reference sequences for {} @SQ lines",
                            reference->NumSequences(), std::size(sequences))};
        }
        for (std::size_t i = 0; i < std::size(sequences); ++i) {
            ReferenceSequence& sq{sequences[i]};
            const auto refId = static_cast<std::int32_t>(i);
            const std::int64_t length{reference->Length(refId)};
            if (length != sq.Length()) {
                throw std::runtime_error{
                    std::format("CramWriter: reference {} has {} bases, its @SQ line says {}",
                                sq.Name(), length, sq.Length())};
            }
            if (!sq.GetTag("M5")) {
                sq.SetTag("M5", CramReferenceMd5(*reference, refId));
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

        PrepareReferences();

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

    // the CRAI entries of a slice, without their file offsets
    std::vector<CraiEntry> SliceCraiEntries(std::span<const BamRecord> sliceRecords,
                                            const CramSliceHeader& sliceHeader) const
    {
        std::vector<CraiEntry> entries;
        if (!config.WriteCrai) {
            return entries;
        }

        const auto containerRefId = sliceHeader.RefSeqId;
        if (containerRefId >= 0) {
            const std::int64_t span = std::max<std::int64_t>(sliceHeader.AlignmentSpan, 1);
            entries.push_back(CraiEntry{.SequenceId = containerRefId,
                                        .AlignmentStart = sliceHeader.AlignmentStart,
                                        .AlignmentSpan = span});
            return entries;
        }

        if (containerRefId == -1) {
            entries.push_back(CraiEntry{.SequenceId = -1, .AlignmentStart = 0, .AlignmentSpan = 1});
            return entries;
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
                entries.push_back(
                    CraiEntry{.SequenceId = -1, .AlignmentStart = 0, .AlignmentSpan = 1});
                continue;
            }

            const std::int64_t alignmentStart = span.MinPos + 1;
            const std::int64_t alignmentSpan = std::max<std::int64_t>(span.MaxEnd - span.MinPos, 1);
            entries.push_back(CraiEntry{.SequenceId = span.SequenceId,
                                        .AlignmentStart = alignmentStart,
                                        .AlignmentSpan = alignmentSpan});
        }
        return entries;
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

    ContainerEncodingPlan BuildContainerEncodingPlan(
        std::span<const std::span<const BamRecord>> slices) const
    {
        ContainerEncodingPlan plan;
        plan.SliceData.resize(std::size(slices));

        plan.CompressionHeader.PreservationMap.ReadNamesIncluded = true;
        plan.CompressionHeader.PreservationMap.ApDelta = false;
        plan.CompressionHeader.PreservationMap.ReferenceRequired = ReferenceBased();

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

        for (std::size_t sliceIndex = 0; sliceIndex < std::size(slices); ++sliceIndex) {
            const auto& sliceRecords = slices[sliceIndex];
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

    static ContainerStats ComputeContainerStats(std::span<const std::span<const BamRecord>> slices)
    {
        if (std::empty(slices)) {
            return {};
        }

        ContainerStats stats;
        std::int32_t minPos = std::numeric_limits<std::int32_t>::max();
        std::int32_t maxEnd = 0;
        bool haveFirstRef = false;
        std::int32_t firstRefId = -2;
        bool allSameRef = true;

        for (const auto& sliceRecords : slices) {
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
        bool allowBlockParallel) const
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

        std::string referenceBuffer;  // the reference bases of a record or of the slice
        for (std::size_t i = 0; i < std::size(records); ++i) {
            EncodeRecord(records[i], recordTagData[i], extStore, containerRefId, referenceBuffer);
        }

        CramSliceHeader sliceHeader;
        sliceHeader.RefSeqId = containerRefId;
        sliceHeader.AlignmentStart = minPos + 1;
        sliceHeader.AlignmentSpan = maxEnd - minPos;
        sliceHeader.NumRecords = numRecords;
        sliceHeader.RecordCounter = recordCounterStart;
        sliceHeader.EmbeddedRefBlockId = -1;
        // The MD5 of the reference bases of a slice on one reference (zero where the slice
        // reaches past the end of the reference: no check).
        if (ReferenceBased() && (containerRefId >= 0) && (maxEnd > minPos)) {
            if (maxEnd <= reference->Length(containerRefId)) {
                sliceHeader.RefMd5 = ReferenceMd5(*reference, containerRefId, minPos,
                                                  maxEnd - minPos, referenceBuffer);
            }
        }

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

    // Encodes slices into one container. The record counters start at 0; WriteEncoded moves
    // them. Without the workers (any thread may call it), with thread-local compressors.
    CramEncodedContainer EncodeSlices(std::span<const std::span<const BamRecord>> sliceRecords,
                                      const bool useWorkers) const
    {
        const ContainerStats containerStats = ComputeContainerStats(sliceRecords);
        const ContainerEncodingPlan encodingPlan = BuildContainerEncodingPlan(sliceRecords);

        std::vector<std::int64_t> sliceRecordCounters(std::size(sliceRecords), 0);
        std::int64_t runningRecordCounter = 0;
        for (std::size_t i = 0; i < std::size(sliceRecords); ++i) {
            sliceRecordCounters[i] = runningRecordCounter;
            runningRecordCounter += static_cast<std::int64_t>(std::size(sliceRecords[i]));
        }

        std::vector<CramSlice> slices(std::size(sliceRecords));
        const bool useParallelSlices = useWorkers && compressionPool &&
                                       config.CompressionWorkers > 1 && std::size(sliceRecords) > 1;
        if (useParallelSlices) {
            std::atomic<std::size_t> nextSlice{0};
            const std::size_t sliceCount = std::size(sliceRecords);
            const auto taskCount = static_cast<std::int32_t>(
                std::ranges::min(config.CompressionWorkers, std::size(sliceRecords)));
            Parallel::Dispatch(
                compressionPool,
                MakeWorkStealingTask(
                    &nextSlice, sliceCount,
                    [this, &slices, &sliceRecords, &encodingPlan, &containerStats,
                     &sliceRecordCounters](std::size_t sliceIndex) {
                        slices[sliceIndex] = EncodeSlice(
                            sliceRecords[sliceIndex], encodingPlan.SliceData[sliceIndex].RecordTags,
                            encodingPlan.DataSeriesMethods, containerStats.RefSeqId,
                            sliceRecordCounters[sliceIndex],
                            /*useThreadLocalGzip=*/true, /*allowBlockParallel=*/false);
                    }),
                taskCount);
        } else {
            for (std::size_t sliceIndex = 0; sliceIndex < std::size(sliceRecords); ++sliceIndex) {
                slices[sliceIndex] = EncodeSlice(
                    sliceRecords[sliceIndex], encodingPlan.SliceData[sliceIndex].RecordTags,
                    encodingPlan.DataSeriesMethods, containerStats.RefSeqId,
                    sliceRecordCounters[sliceIndex],
                    /*useThreadLocalGzip=*/!useWorkers, /*allowBlockParallel=*/useWorkers);
            }
        }

        CramEncodedContainer encoded;
        CramContainer& container{encoded.Container};
        container.Header.RefSeqId = containerStats.RefSeqId;
        container.Header.StartPos = containerStats.StartPos;
        container.Header.AlignmentSpan = containerStats.AlignmentSpan;
        container.Header.NumRecords = containerStats.NumRecords;
        container.Header.RecordCounter = 0;
        container.Header.Bases = containerStats.Bases;
        container.CompressionHeader = encodingPlan.CompressionHeader;
        container.Slices = std::move(slices);
        encoded.SliceIndex.reserve(std::size(sliceRecords));
        for (std::size_t i = 0; i < std::size(sliceRecords); ++i) {
            encoded.SliceIndex.push_back(
                SliceCraiEntries(sliceRecords[i], container.Slices[i].Header));
        }
        return encoded;
    }

    // Writes an encoded container after the containers written before: the record counters
    // continue, and the CRAI entries get their file offsets.
    void WriteEncoded(CramEncodedContainer encoded)
    {
        CramContainer& container{encoded.Container};
        if (std::size(encoded.SliceIndex) != std::size(container.Slices)) {
            throw std::runtime_error("CramWriter: encoded container without its slice index");
        }
        container.Header.RecordCounter = globalRecordCounter;
        for (CramSlice& slice : container.Slices) {
            slice.Header.RecordCounter += globalRecordCounter;
        }

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
            for (const CraiEntry& entry : encoded.SliceIndex[i]) {
                AddCraiEntry(entry.SequenceId, entry.AlignmentStart, entry.AlignmentSpan,
                             containerOffset, sliceOffset, sliceSize);
            }
        }

        QueueContainerWrite(std::move(serialized));
        globalRecordCounter += container.Header.NumRecords;
    }

    void FlushContainer()
    {
        if (std::empty(pendingSlices)) {
            return;
        }
        const std::vector<std::span<const BamRecord>> slices(std::begin(pendingSlices),
                                                             std::end(pendingSlices));
        WriteEncoded(EncodeSlices(slices, /*useWorkers=*/true));
        pendingSlices.clear();
    }

    // the records written one by one go before a container of EncodeContainer
    void FlushPendingRecords()
    {
        if (!std::empty(pendingRecords)) {
            pendingSlices.push_back(std::move(pendingRecords));
            pendingRecords.clear();
            pendingRecords.reserve(static_cast<std::size_t>(config.RecordsPerSlice));
        }
        FlushContainer();
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
                      CramExternalBlockStore& extStore, std::int32_t containerRefId,
                      std::string& referenceBuffer) const
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
            // The features go to their own external blocks, so FN can follow them.
            // The reference bases under the alignment (fewer past the end of the reference),
            // from windowStart.
            const bool againstReference{HasReference(record)};
            std::string_view window;
            std::int64_t windowStart{0};
            if (againstReference) {
                const std::int64_t length{reference->Length(record.RefId())};
                windowStart = std::clamp<std::int64_t>(record.Pos(), 0, length);
                const std::int64_t windowEnd{std::clamp<std::int64_t>(
                    record.Pos() + ReferenceSpan(record), windowStart, length)};
                window = reference->Bases(record.RefId(), windowStart, windowEnd - windowStart,
                                          referenceBuffer);
            }
            std::int32_t numFeatures = 0;
            std::int32_t readPos = 1;
            std::int32_t seqOffset = 0;
            std::int64_t refPos = record.Pos();
            std::int32_t basesCovered = 0;
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
                            basesCovered += static_cast<std::int32_t>(std::size(data));
                            if (againstReference) {
                                WriteBasesAgainstReference(extStore, data, window,
                                                           refPos - windowStart, readPos,
                                                           previousFeaturePos, numFeatures);
                            } else {
                                WriteFeatureHeader(extStore, 'b', readPos, previousFeaturePos);
                                WriteByteArrayLen(extStore, BB_BLOCK_ID, data);
                                ++numFeatures;
                                wroteFeature = true;
                            }
                        }
                        readPos += length;
                        seqOffset += length;
                        refPos += length;
                        break;
                    }
                    case CigarOpType::I: {
                        const auto data = QuerySpan(seq, seqOffset, length);
                        if (!std::empty(data)) {
                            basesCovered += static_cast<std::int32_t>(std::size(data));
                            WriteFeatureHeader(extStore, 'I', readPos, previousFeaturePos);
                            WriteByteArrayStop(extStore, IN_BLOCK_ID, data, std::byte{0});
                            ++numFeatures;
                            wroteFeature = true;
                        }
                        readPos += length;
                        seqOffset += length;
                        break;
                    }
                    case CigarOpType::S: {
                        const auto data = QuerySpan(seq, seqOffset, length);
                        if (!std::empty(data)) {
                            basesCovered += static_cast<std::int32_t>(std::size(data));
                            WriteFeatureHeader(extStore, 'S', readPos, previousFeaturePos);
                            WriteByteArrayStop(extStore, SC_BLOCK_ID, data, std::byte{0});
                            ++numFeatures;
                            wroteFeature = true;
                        }
                        readPos += length;
                        seqOffset += length;
                        break;
                    }
                    case CigarOpType::D:
                        WriteFeatureHeader(extStore, 'D', readPos, previousFeaturePos);
                        extStore.WriteItf8(DL_BLOCK_ID, length);
                        ++numFeatures;
                        wroteFeature = true;
                        refPos += length;
                        break;
                    case CigarOpType::N:
                        WriteFeatureHeader(extStore, 'N', readPos, previousFeaturePos);
                        extStore.WriteItf8(RS_BLOCK_ID, length);
                        ++numFeatures;
                        wroteFeature = true;
                        refPos += length;
                        break;
                    case CigarOpType::H:
                        WriteFeatureHeader(extStore, 'H', readPos, previousFeaturePos);
                        extStore.WriteItf8(HC_BLOCK_ID, length);
                        ++numFeatures;
                        wroteFeature = true;
                        break;
                    case CigarOpType::P:
                        WriteFeatureHeader(extStore, 'P', readPos, previousFeaturePos);
                        extStore.WriteItf8(PD_BLOCK_ID, length);
                        ++numFeatures;
                        wroteFeature = true;
                        break;
                }
            }

            // A CIGAR without the bases of the read: store them as one run.
            if (!wroteFeature && (basesCovered == 0) && !std::empty(seq)) {
                WriteFeatureHeader(extStore, 'b', 1, previousFeaturePos);
                const auto* seqPtr = reinterpret_cast<const std::byte*>(seq.data());
                WriteByteArrayLen(extStore, BB_BLOCK_ID, {seqPtr, std::size(seq)});
                ++numFeatures;
            }
            extStore.WriteItf8(FN_BLOCK_ID, numFeatures);

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

CramEncodedContainer CramWriter::EncodeContainer(const std::span<const BamRecord> records) const
{
    std::vector<std::span<const BamRecord>> slices;
    const auto perSlice = static_cast<std::size_t>(impl_->config.RecordsPerSlice);
    for (std::size_t first = 0; first < std::size(records); first += perSlice) {
        slices.push_back(records.subspan(first, std::min(perSlice, std::size(records) - first)));
    }
    return impl_->EncodeSlices(slices, /*useWorkers=*/false);
}

void CramWriter::WriteContainer(CramEncodedContainer container)
{
    if (std::empty(container.Container.Slices)) {
        return;
    }
    impl_->FlushPendingRecords();
    impl_->WriteEncoded(std::move(container));
}

std::string CramReferenceMd5(const std::string_view sequence)
{
    return detail::Md5DigestToHex(detail::ComputeUppercaseMd5(sequence));
}

std::string CramReferenceMd5(const CramReference& reference, const std::int32_t refId)
{
    std::string buffer;
    return detail::Md5DigestToHex(
        ReferenceMd5(reference, refId, 0, reference.Length(refId), buffer));
}

}  // namespace Samoa
}  // namespace PacBio
