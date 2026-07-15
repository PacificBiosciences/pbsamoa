#include <pbsamoa/io/CramReader.hpp>

#include <pbsamoa/cram/CramCodec.hpp>
#include <pbsamoa/cram/CramCompression.hpp>
#include <pbsamoa/cram/CramStructs.hpp>

#include "BinaryUtils.hpp"
#include "CramInternal.hpp"
#include "CramStructs.hpp"
#include "Md5.hpp"
#include "ReaderUtils.hpp"

#include <pbcopper/parallel/ThreadPool.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <format>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include <cctype>
#include <cstdlib>

namespace PacBio {
namespace Samoa {

namespace {

struct DecodedFeature
{
    char Code{};
    std::int32_t ReadPos{0};  // absolute, 1-based
    std::int32_t Length{0};
    bool IsSubstitutionCode{false};
    std::vector<std::byte> Data;
    std::optional<std::uint8_t> Quality;  // 'B': QS byte placed at ReadPos-1
};

std::size_t Itf8EncodedLength(std::uint8_t b0)
{
    if ((b0 & 0x80) == 0) {
        return 1;
    }
    if ((b0 & 0xC0) == 0x80) {
        return 2;
    }
    if ((b0 & 0xE0) == 0xC0) {
        return 3;
    }
    if ((b0 & 0xF0) == 0xE0) {
        return 4;
    }
    return 5;
}

std::size_t Ltf8EncodedLength(std::uint8_t b0)
{
    if ((b0 & 0x80) == 0) {
        return 1;
    }
    if ((b0 & 0xC0) == 0x80) {
        return 2;
    }
    if ((b0 & 0xE0) == 0xC0) {
        return 3;
    }
    if ((b0 & 0xF0) == 0xE0) {
        return 4;
    }
    if ((b0 & 0xF8) == 0xF0) {
        return 5;
    }
    if ((b0 & 0xFC) == 0xF8) {
        return 6;
    }
    if ((b0 & 0xFE) == 0xFC) {
        return 7;
    }
    if (b0 == 0xFE) {
        return 8;
    }
    return 9;
}

std::vector<std::vector<TagTriple>> ParseTagDictionary(std::span<const std::byte> dictionaryBytes)
{
    std::vector<std::vector<TagTriple>> dictionary;
    std::vector<TagTriple> current;

    for (std::size_t i = 0; i < std::size(dictionaryBytes);) {
        if (dictionaryBytes[i] == std::byte{0}) {
            dictionary.push_back(std::move(current));
            current.clear();
            ++i;
            continue;
        }
        if ((i + 2) >= std::size(dictionaryBytes)) {
            break;
        }
        current.push_back(TagTriple{
            static_cast<char>(static_cast<std::uint8_t>(dictionaryBytes[i])),
            static_cast<char>(static_cast<std::uint8_t>(dictionaryBytes[i + 1])),
            static_cast<char>(static_cast<std::uint8_t>(dictionaryBytes[i + 2])),
        });
        i += 3;
    }

    if (!std::empty(current)) {
        dictionary.push_back(std::move(current));
    }
    if (std::empty(dictionary)) {
        dictionary.emplace_back();  // index 0: empty tag set
    }

    return dictionary;
}

char ByteToChar(std::byte value) { return static_cast<char>(std::to_integer<std::uint8_t>(value)); }

std::uint8_t ByteToUInt8(std::byte value) { return std::to_integer<std::uint8_t>(value); }

bool IsZeroByte(std::byte value) { return value == std::byte{0}; }

bool IsNonZeroByte(std::byte value) { return value != std::byte{0}; }

template <typename T, typename Converter>
void CopyBytesTo(std::span<const std::byte> source, std::span<T> destination,
                 std::size_t destOffset, Converter convert)
{
    const std::size_t sourceSize{std::size(source)};
    const std::size_t destinationSize{std::size(destination)};
    if (destOffset >= destinationSize) {
        return;
    }

    const std::size_t copyLen{std::ranges::min(sourceSize, destinationSize - destOffset)};
    for (std::size_t i{0}; i < copyLen; ++i) {
        destination[destOffset + i] = convert(source[i]);
    }
}

void CopyBytesToString(std::span<const std::byte> source, std::string& destination,
                       std::size_t destOffset)
{
    CopyBytesTo(source, std::span<char>{destination}, destOffset, ByteToChar);
}

void CopyBytesToVector(std::span<const std::byte> source, std::vector<std::uint8_t>& destination,
                       std::size_t destOffset)
{
    CopyBytesTo(source, std::span<std::uint8_t>{destination}, destOffset, ByteToUInt8);
}

libdeflate_decompressor* ThreadLocalGzipDecompressor()
{
    const thread_local CramGzipDecompressorContext THREAD_LOCAL_GZIP_DECOMPRESSOR_CONTEXT{};
    if (!THREAD_LOCAL_GZIP_DECOMPRESSOR_CONTEXT.Decompressor) {
        throw std::runtime_error("CramReader: failed to allocate gzip decompressor");
    }
    return THREAD_LOCAL_GZIP_DECOMPRESSOR_CONTEXT.Decompressor.get();
}

void AppendCigarOp(std::vector<CigarOp>& cigar, CigarOpType type, std::int32_t length)
{
    if (length <= 0) {
        return;
    }
    const std::uint32_t ulen = length;
    if (!std::empty(cigar) && cigar.back().Type() == type) {
        const auto mergedLen = cigar.back().Length() + ulen;
        cigar.back() = CigarOp{type, mergedLen};
        return;
    }
    cigar.push_back(CigarOp{type, ulen});
}

RawRecord ToRawRecord(const BamRecord& record)
{
    const std::vector<std::byte> bytes = record.SerializeToBam();
    return RawRecord{std::span<const std::byte>{bytes}};
}

template <typename T>
void AppendLittleEndian(std::vector<std::byte>& payload, T value)
{
    const auto* ptr = reinterpret_cast<const std::byte*>(&value);
    payload.insert(std::end(payload), ptr, ptr + sizeof(T));
}

std::vector<std::byte> RequirePayloadWidth(char type, std::vector<std::byte> payload,
                                           std::size_t expectedWidth)
{
    if (std::size(payload) != expectedWidth) {
        throw std::runtime_error{
            std::format("CramReader: tag type '{}' requires {} payload byte(s), got {}", type,
                        expectedWidth, std::size(payload))};
    }
    return payload;
}

std::vector<std::byte> DecodeRequiredByteArrayPayload(char type, CramCodec& codec,
                                                      CramBitReader& coreReader,
                                                      CramExternalBlockStore& extStore)
{
    if (codec.DecodeKind() != CramCodecDecodeKind::BYTE_ARRAY) {
        throw std::runtime_error{
            std::format("CramReader: tag type '{}' requires byte-array codec decode", type)};
    }
    return codec.DecodeByteArray(coreReader, extStore);
}

std::vector<std::byte> DecodeRequiredFixedWidthPayload(char type, CramCodec& codec,
                                                       CramBitReader& coreReader,
                                                       CramExternalBlockStore& extStore,
                                                       std::size_t expectedWidth)
{
    return RequirePayloadWidth(
        type, DecodeRequiredByteArrayPayload(type, codec, coreReader, extStore), expectedWidth);
}

void DecodeReadName(std::string& readName, CramCodec& codec, CramBitReader& coreReader,
                    CramExternalBlockStore& extStore)
{
    const auto nameView = codec.DecodeByteArrayView(coreReader, extStore);
    readName.assign(reinterpret_cast<const char*>(nameView.data()), std::size(nameView));
}

void DecodeQualityArray(std::vector<std::uint8_t>& qualities, CramCodec& codec,
                        CramBitReader& coreReader, CramExternalBlockStore& extStore,
                        std::int32_t readLength)
{
    qualities.reserve(readLength);
    for (std::int32_t qi = 0; qi < readLength; ++qi) {
        qualities.push_back(static_cast<std::uint8_t>(codec.DecodeByte(coreReader, extStore)));
    }
}

std::int32_t FeatureLengthOrOne(const DecodedFeature& feature)
{
    if (feature.Length > 0) {
        return feature.Length;
    }
    return 1;
}

void AppendBytes(std::vector<std::byte>& destination, std::span<const std::byte> source)
{
    destination.insert(std::end(destination), std::begin(source), std::end(source));
}

std::int32_t DecodeEncodedItf8(std::span<const std::byte> encoded)
{
    std::size_t bytesRead = 0;
    return ReadItf8(encoded, bytesRead);
}

void DecodeSubstitutionMatrixRow(std::array<std::array<char, 4>, 5>& matrix, std::size_t row,
                                 std::array<char, 4> symbols, std::uint8_t packed)
{
    for (int i = 0; i < 4; ++i) {
        const std::size_t idx = (packed >> (6 - 2 * i)) & 0x03;
        matrix[row][idx] = symbols[static_cast<std::size_t>(i)];
    }
}

// Handles tag type 'A' (char): BYTE_ARRAY codec path uses DecodeByte, not DecodeInt.
std::vector<std::byte> DecodeByteTagPayload(char type, CramCodec& codec, CramBitReader& coreReader,
                                            CramExternalBlockStore& extStore)
{
    if (codec.DecodeKind() == CramCodecDecodeKind::BYTE_ARRAY) {
        return DecodeRequiredFixedWidthPayload(type, codec, coreReader, extStore, 1);
    }
    std::vector<std::byte> payload;
    payload.push_back(codec.DecodeByte(coreReader, extStore));
    return payload;
}

// Handles the 7 numeric tag types (c C s S i I f): scalar path uses DecodeInt.
template <typename T>
std::vector<std::byte> DecodeIntTagPayload(char type, std::size_t width, CramCodec& codec,
                                           CramBitReader& coreReader,
                                           CramExternalBlockStore& extStore)
{
    if (codec.DecodeKind() == CramCodecDecodeKind::BYTE_ARRAY) {
        return DecodeRequiredFixedWidthPayload(type, codec, coreReader, extStore, width);
    }
    std::vector<std::byte> payload;
    const T value = codec.DecodeInt(coreReader, extStore);
    AppendLittleEndian(payload, value);
    return payload;
}

std::vector<std::byte> DecodeTagPayload(char type, CramCodec& codec, CramBitReader& coreReader,
                                        CramExternalBlockStore& extStore)
{
    const auto decodeKind = codec.DecodeKind();
    std::vector<std::byte> payload;

    switch (type) {
        case 'A':
            return DecodeByteTagPayload(type, codec, coreReader, extStore);
        case 'c':
            return DecodeIntTagPayload<std::int8_t>(type, 1, codec, coreReader, extStore);
        case 'C':
            return DecodeIntTagPayload<std::uint8_t>(type, 1, codec, coreReader, extStore);
        case 's':
            return DecodeIntTagPayload<std::int16_t>(type, 2, codec, coreReader, extStore);
        case 'S':
            return DecodeIntTagPayload<std::uint16_t>(type, 2, codec, coreReader, extStore);
        case 'i':
            return DecodeIntTagPayload<std::int32_t>(type, 4, codec, coreReader, extStore);
        case 'I':
            return DecodeIntTagPayload<std::uint32_t>(type, 4, codec, coreReader, extStore);
        case 'f':
            return DecodeIntTagPayload<std::uint32_t>(type, sizeof(float), codec, coreReader,
                                                      extStore);
        case 'Z':
        case 'H':
            payload = DecodeRequiredByteArrayPayload(type, codec, coreReader, extStore);
            if (std::empty(payload) || (payload.back() != std::byte{0})) {
                payload.push_back(std::byte{0});
            }
            return payload;
        case 'B':
            return DecodeRequiredByteArrayPayload(type, codec, coreReader, extStore);
        default:
            if (decodeKind != CramCodecDecodeKind::BYTE_ARRAY) {
                throw std::runtime_error{
                    std::format("CramReader: unsupported scalar tag type '{}'", type)};
            }
            return DecodeRequiredByteArrayPayload(type, codec, coreReader, extStore);
    }
}

template <typename Key>
CramCodec* LookupSortedCodec(std::span<const std::pair<Key, CramCodec*>> codecs, const Key& key)
{
    const auto it = std::ranges::lower_bound(codecs, key, {}, &std::pair<Key, CramCodec*>::first);
    if (it != std::end(codecs) && it->first == key) {
        return it->second;
    }
    return nullptr;
}

std::array<std::array<char, 4>, 5> BuildSubstitutionMatrixLookup(
    std::span<const std::byte, 5> encoded)
{
    std::array<std::array<char, 4>, 5> matrix{{
        {{'C', 'G', 'T', 'N'}},  // A
        {{'A', 'G', 'T', 'N'}},  // C
        {{'A', 'C', 'T', 'N'}},  // G
        {{'A', 'C', 'G', 'N'}},  // T
        {{'A', 'C', 'G', 'T'}},  // N
    }};

    if (std::ranges::all_of(encoded, IsZeroByte)) {
        return matrix;
    }

    DecodeSubstitutionMatrixRow(matrix, 0, {'C', 'G', 'T', 'N'}, ByteToUInt8(encoded[0]));
    DecodeSubstitutionMatrixRow(matrix, 1, {'A', 'G', 'T', 'N'}, ByteToUInt8(encoded[1]));
    DecodeSubstitutionMatrixRow(matrix, 2, {'A', 'C', 'T', 'N'}, ByteToUInt8(encoded[2]));
    DecodeSubstitutionMatrixRow(matrix, 3, {'A', 'C', 'G', 'N'}, ByteToUInt8(encoded[3]));
    DecodeSubstitutionMatrixRow(matrix, 4, {'A', 'C', 'G', 'T'}, ByteToUInt8(encoded[4]));

    return matrix;
}

std::size_t RefBaseRowIndex(char base)
{
    switch (base) {
        case 'A':
        case 'a':
            return 0;
        case 'C':
        case 'c':
            return 1;
        case 'G':
        case 'g':
            return 2;
        case 'T':
        case 't':
            return 3;
        default:
            return 4;
    }
}

bool HasNonZeroRefMd5(std::span<const std::byte, 16> refMd5)
{
    return std::ranges::any_of(refMd5, IsNonZeroByte);
}

bool SliceRequiresReferenceMd5Validation(const CramSliceHeader& sliceHeader, bool referenceRequired)
{
    return referenceRequired && sliceHeader.RefSeqId >= 0 && HasNonZeroRefMd5(sliceHeader.RefMd5);
}

std::string SliceReferenceLabel(const CramSliceHeader& sliceHeader, const SamHeader& header)
{
    const auto refId = sliceHeader.RefSeqId;
    if (refId >= 0 && static_cast<std::size_t>(refId) < std::size(header.ReferenceSequences())) {
        const auto& reference = header.ReferenceSequences()[static_cast<std::size_t>(refId)];
        return std::format("{} (id={})", std::string{reference.Name()}, refId);
    }
    return std::format("id={}", refId);
}

std::string SliceSpanLabel(const CramSliceHeader& sliceHeader)
{
    const std::int64_t start = sliceHeader.AlignmentStart;
    const std::int64_t span = sliceHeader.AlignmentSpan;
    if (start > 0 && span > 0) {
        return std::format("{}-{}", start, start + span - 1);
    }
    return std::format("start={} span={}", start, span);
}

std::vector<std::byte> ResolveSliceReferenceWindow(
    const CramSliceHeader& sliceHeader, std::span<const std::byte> embeddedReference,
    const std::vector<std::string>& externalReferenceById, const SamHeader& header)
{
    const auto referenceLabel = SliceReferenceLabel(sliceHeader, header);
    const auto spanLabel = SliceSpanLabel(sliceHeader);

    if (sliceHeader.AlignmentStart <= 0 || sliceHeader.AlignmentSpan <= 0) {
        throw std::runtime_error{
            std::format("CramReader: reference MD5 validation requires positive "
                        "slice coordinates for "
                        "reference {} span {}",
                        referenceLabel, spanLabel)};
    }

    const std::size_t windowLength = sliceHeader.AlignmentSpan;
    const std::size_t windowStart = sliceHeader.AlignmentStart - 1;

    if (!std::empty(embeddedReference)) {
        if (std::size(embeddedReference) < windowLength) {
            throw std::runtime_error{
                std::format("CramReader: reference MD5 validation required but "
                            "embedded reference is too "
                            "short for reference {} span {} (have {}, need {})",
                            referenceLabel, spanLabel, std::size(embeddedReference), windowLength)};
        }
        const auto window = embeddedReference.first(windowLength);
        return {std::begin(window), std::end(window)};
    }

    if (sliceHeader.RefSeqId >= 0 &&
        static_cast<std::size_t>(sliceHeader.RefSeqId) < std::size(externalReferenceById)) {
        const auto& fullReference =
            externalReferenceById[static_cast<std::size_t>(sliceHeader.RefSeqId)];
        if (!fullReference.empty()) {
            if (windowStart >= std::size(fullReference) ||
                windowLength > (std::size(fullReference) - windowStart)) {
                throw std::runtime_error{
                    std::format("CramReader: reference MD5 validation required but "
                                "FASTA reference is too "
                                "short for reference {} span {} (length {})",
                                referenceLabel, spanLabel, std::size(fullReference))};
            }

            const auto src = std::string_view{fullReference}.substr(windowStart, windowLength);
            std::vector<std::byte> window(windowLength);
            std::ranges::transform(src, window.begin(), [](char c) {
                return static_cast<std::byte>(static_cast<unsigned char>(c));
            });
            return window;
        }
    }

    throw std::runtime_error{
        std::format("CramReader: reference MD5 validation required but no "
                    "reference bases are available for "
                    "reference {} span {}",
                    referenceLabel, spanLabel)};
}

std::vector<std::byte> UppercaseBases(std::span<const std::byte> bases)
{
    std::vector<std::byte> normalized(std::size(bases));
    std::ranges::transform(bases, normalized.begin(), [](std::byte b) {
        return static_cast<std::byte>(std::toupper(std::to_integer<unsigned char>(b)));
    });
    return normalized;
}

void ValidateSliceReferenceMd5(const CramSliceHeader& sliceHeader, bool referenceRequired,
                               std::span<const std::byte> embeddedReference,
                               const std::vector<std::string>& externalReferenceById,
                               const SamHeader& header)
{
    if (!SliceRequiresReferenceMd5Validation(sliceHeader, referenceRequired)) {
        return;
    }

    const auto referenceLabel = SliceReferenceLabel(sliceHeader, header);
    const auto spanLabel = SliceSpanLabel(sliceHeader);

    const auto referenceWindow =
        ResolveSliceReferenceWindow(sliceHeader, embeddedReference, externalReferenceById, header);
    const auto normalizedReference = UppercaseBases(referenceWindow);
    const auto observedMd5 = detail::ComputeMd5(normalizedReference);
    if (observedMd5 != sliceHeader.RefMd5) {
        throw std::runtime_error{
            std::format("CramReader: reference MD5 mismatch for reference {} span "
                        "{} (expected {}, observed "
                        "{})",
                        referenceLabel, spanLabel, detail::Md5DigestToHex(sliceHeader.RefMd5),
                        detail::Md5DigestToHex(observedMd5))};
    }
}

std::unordered_map<std::string, std::string> LoadFastaSequences(const std::filesystem::path& path)
{
    std::ifstream in{path};
    if (!in) {
        throw std::runtime_error{
            std::format("CramReader: cannot open reference FASTA {}", path.string())};
    }

    std::unordered_map<std::string, std::string> references;
    std::string currentName;
    std::string line;
    while (std::getline(in, line)) {
        if (!std::empty(line) && line.back() == '\r') {
            line.pop_back();
        }
        if (std::empty(line)) {
            continue;
        }
        if (line[0] == '>') {
            line.erase(0, 1);
            const auto wsPos = line.find_first_of(" \t");
            currentName = line.substr(0, wsPos);
            references.try_emplace(currentName);
            continue;
        }
        if (std::empty(currentName)) {
            continue;
        }

        auto& seq = references[currentName];
        for (const unsigned char c : line) {
            if (std::isspace(c) == 0) {
                seq.push_back(static_cast<char>(std::toupper(c)));
            }
        }
    }
    return references;
}

void ResolveDownstreamMates(std::vector<BamRecord>& records,
                            std::span<const std::int32_t> downstreamMateDistance)
{
    const auto n = std::size(records);
    if (std::size(downstreamMateDistance) != n) {
        return;
    }

    for (std::size_t i = 0; i < n; ++i) {
        const auto dist = downstreamMateDistance[i];
        if (dist <= 0) {
            continue;
        }
        const auto mateIndex = i + static_cast<std::size_t>(dist);
        if (mateIndex >= n) {
            continue;
        }

        auto& rec = records[i];
        auto& mate = records[mateIndex];

        rec.NextRefId(mate.RefId()).NextPos(mate.Pos());
        mate.NextRefId(rec.RefId()).NextPos(rec.Pos());

        if (rec.RefId() < 0 || mate.RefId() < 0 || rec.RefId() != mate.RefId()) {
            rec.Tlen(0);
            mate.Tlen(0);
            continue;
        }

        const auto leftPos = std::ranges::min(rec.Pos(), mate.Pos());
        const auto rightEnd = std::ranges::max(rec.ReferenceEnd(), mate.ReferenceEnd());
        const auto templateLen = rightEnd - leftPos;
        if (templateLen <= 0) {
            continue;
        }

        if (rec.Pos() <= mate.Pos()) {
            rec.Tlen(templateLen);
            mate.Tlen(-templateLen);
        } else {
            rec.Tlen(-templateLen);
            mate.Tlen(templateLen);
        }
    }
}

struct SliceQueryCandidate
{
    std::int64_t ContainerOffset{};
    std::int64_t SliceOffset{};
    std::int64_t SliceSize{};
};

std::array<std::int64_t, 3> SliceQueryCandidateKey(const SliceQueryCandidate& candidate)
{
    return {candidate.ContainerOffset, candidate.SliceOffset, candidate.SliceSize};
}

std::vector<SliceQueryCandidate> BuildSliceCandidates(const CraiIndex& index, std::int32_t refId,
                                                      std::int32_t beg, std::int32_t end)
{
    std::vector<SliceQueryCandidate> candidates;

    const auto entries = index.EntriesForReference(refId);
    candidates.reserve(std::size(entries));

    for (const CraiEntry& entry : entries) {
        if (refId >= 0) {
            const std::int64_t rowBeg = entry.AlignmentStart - 1;
            const std::int64_t rowEnd =
                rowBeg + std::max<std::int64_t>(entry.AlignmentSpan, static_cast<std::int64_t>(1));
            if (!(rowBeg < end && rowEnd > beg)) {
                continue;
            }
        } else if (entry.SequenceId != -1) {
            continue;
        }

        candidates.push_back(SliceQueryCandidate{
            .ContainerOffset = entry.ContainerOffset,
            .SliceOffset = entry.SliceOffset,
            .SliceSize = entry.SliceSize,
        });
    }

    std::ranges::sort(candidates, {}, SliceQueryCandidateKey);
    const auto duplicateStart =
        std::ranges::unique(candidates, std::ranges::equal_to{}, SliceQueryCandidateKey);
    candidates.erase(duplicateStart.begin(), std::end(candidates));
    return candidates;
}

std::int64_t QueryRecordEnd(const BamRecord& record)
{
    const std::int64_t pos{record.Pos()};
    return NonEmptyAlignmentEnd(pos, static_cast<std::int64_t>(record.ReferenceEnd()));
}

bool KeepRecordForQuery(const BamRecord& record, std::int32_t refId, std::int32_t beg,
                        std::int32_t end)
{
    if (refId == -1) {
        return record.RefId() < 0;
    }

    if (record.RefId() != refId) {
        return false;
    }
    return record.Pos() < end && QueryRecordEnd(record) > beg;
}

}  // namespace

// ---------------------------------------------------------------------------
// Impl
// ---------------------------------------------------------------------------

struct CramReader::Impl
{
    struct DecodeRecordContext
    {
        CramBitReader& CoreReader;
        CramExternalBlockStore& ExtStore;
        std::span<const std::pair<CramDataSeries, CramCodec*>> Codecs;
        std::span<const std::pair<std::int32_t, CramCodec*>> TagCodecs;
        const std::vector<std::vector<TagTriple>>& TagDictionary;
        const CramSliceHeader& SliceHeader;
        const std::array<std::array<char, 4>, 5>& SubstitutionMatrix;
        std::span<const std::byte> EmbeddedReference;
        std::int64_t& PreviousApos;
        std::int32_t& DownstreamMateDistance;
        bool ApDelta{false};
        bool ReadNamesIncluded{false};
        const std::vector<std::string>& ExternalReferenceById;
    };

    std::filesystem::path path;
    CramReaderConfig config;
    std::ifstream file;
    CramFileDefinition fileDef;
    SamHeader header;
    bool atEof{false};

    // Current container state
    std::vector<BamRecord> pendingRecords;
    std::size_t pendingIdx{0};
    std::vector<std::string> externalReferenceById;
    CramGzipDecompressorContext gzipContext;
    std::shared_ptr<Parallel::ThreadPool<>> decompressionPool;

    static std::optional<char> ReferenceBaseAt(const DecodeRecordContext& ctx, std::int32_t refId,
                                               std::int32_t refPos1Based)
    {
        if (!std::empty(ctx.EmbeddedReference) && ctx.SliceHeader.RefSeqId >= 0 && refId >= 0 &&
            refId == ctx.SliceHeader.RefSeqId && ctx.SliceHeader.AlignmentStart > 0 &&
            refPos1Based >= ctx.SliceHeader.AlignmentStart) {
            const auto offset =
                static_cast<std::size_t>(refPos1Based - ctx.SliceHeader.AlignmentStart);
            if (offset < std::size(ctx.EmbeddedReference)) {
                return static_cast<char>(static_cast<std::uint8_t>(ctx.EmbeddedReference[offset]));
            }
        }

        if (refId >= 0 && static_cast<std::size_t>(refId) < std::size(ctx.ExternalReferenceById) &&
            !ctx.ExternalReferenceById[static_cast<std::size_t>(refId)].empty() &&
            refPos1Based > 0) {
            const std::size_t idx = refPos1Based - 1;
            const auto& ref = ctx.ExternalReferenceById[static_cast<std::size_t>(refId)];
            if (idx < std::size(ref)) {
                return ref[idx];
            }
        }
        return std::nullopt;
    }

    static CramCodec* LookupCodec(const DecodeRecordContext& ctx, const CramDataSeries dataSeries)
    {
        return LookupSortedCodec(ctx.Codecs, dataSeries);
    }

    static CramCodec* LookupTagCodec(const DecodeRecordContext& ctx, const std::int32_t contentId)
    {
        return LookupSortedCodec(ctx.TagCodecs, contentId);
    }

    static char SubstitutionBase(const DecodeRecordContext& ctx, std::int32_t refId,
                                 std::uint8_t code, std::int32_t refPos1Based)
    {
        const auto refBase = ReferenceBaseAt(ctx, refId, refPos1Based).value_or('N');
        return ctx.SubstitutionMatrix[RefBaseRowIndex(refBase)][code & 0x03];
    }

    static void FillReferenceMatches(std::string& sequence, const DecodeRecordContext& ctx,
                                     std::int32_t refId, std::int32_t readPos1Based,
                                     std::int32_t refPos1Based, std::int32_t length)
    {
        for (std::int32_t i = 0; i < length; ++i) {
            const std::size_t readIndex = readPos1Based - 1 + i;
            if (readIndex >= std::size(sequence)) {
                break;
            }
            if (const auto refBase = ReferenceBaseAt(ctx, refId, refPos1Based + i); refBase) {
                sequence[readIndex] = *refBase;
            }
        }
    }

    static void WriteFeatureData(std::string& sequence, std::int32_t readPos1Based,
                                 std::span<const std::byte> dataBytes)
    {
        CopyBytesToString(dataBytes, sequence, readPos1Based - 1);
    }

    void Open()
    {
        if (!gzipContext.Decompressor) {
            throw std::runtime_error("CramReader: failed to allocate gzip decompressor");
        }

        file.open(path, std::ios::binary);
        if (!file) {
            throw std::runtime_error{std::format("CramReader: cannot open {}", path.string())};
        }

        // Read file definition (26 bytes)
        std::vector<std::byte> defBuf(26);
        file.read(reinterpret_cast<char*>(defBuf.data()), 26);
        if (!file || file.gcount() != 26) {
            throw std::runtime_error("CramReader: failed to read file definition");
        }
        fileDef = ParseFileDefinition(defBuf);

        if (fileDef.MajorVersion != 3) {
            throw std::runtime_error{std::format("CramReader: unsupported CRAM version {}.{}",
                                                 fileDef.MajorVersion, fileDef.MinorVersion)};
        }

        // Read header container
        ReadHeaderContainer();
        LoadExternalReferences();
        if (config.DecompressionWorkers > 0) {
            decompressionPool =
                std::make_shared<Parallel::ThreadPool<>>(Parallel::ThreadPool<>::Config{
                    .NumThreads = config.DecompressionWorkers,
                });
        }
    }

    void LoadExternalReferences()
    {
        std::filesystem::path referencePath = config.ReferencePath;
        if (referencePath.empty()) {
            if (const char* envRef = std::getenv("PBSAMOA_CRAM_REFERENCE"); envRef) {
                referencePath = envRef;
            }
        }
        if (referencePath.empty()) {
            return;
        }

        auto referencesByName = LoadFastaSequences(referencePath);
        externalReferenceById.clear();
        externalReferenceById.resize(std::size(header.ReferenceSequences()));

        for (std::size_t i = 0; i < std::size(header.ReferenceSequences()); ++i) {
            const auto& sq = header.ReferenceSequences()[i];
            auto it = referencesByName.find(std::string{sq.Name()});
            if (it != std::end(referencesByName)) {
                externalReferenceById[i] = std::move(it->second);
            }
        }
    }

    std::vector<std::byte> ReadBytes(std::size_t n)
    {
        std::vector<std::byte> buf(n);
        file.clear();  // reset stream error flags before reading
        file.read(reinterpret_cast<char*>(buf.data()), static_cast<std::streamsize>(n));
        const std::size_t got = file.gcount();
        if (got < n) {
            buf.resize(got);
        }
        return buf;
    }

    std::vector<std::byte> ReadExactBytes(std::size_t n, std::string_view context)
    {
        const auto data = ReadBytes(n);
        if (std::size(data) != n) {
            throw std::runtime_error{
                std::format("CramReader: truncated {} (expected {} bytes, got {})", context, n,
                            std::size(data))};
        }
        return data;
    }

    template <typename LenFn>
    std::vector<std::byte> ReadVarIntBytes(LenFn lenFn, std::string_view context)
    {
        auto encoded = ReadExactBytes(1, context);
        const auto totalLen = lenFn(static_cast<std::uint8_t>(encoded[0]));
        if (totalLen > 1) {
            auto tail = ReadExactBytes(totalLen - 1, context);
            encoded.insert(std::end(encoded), std::begin(tail), std::end(tail));
        }
        return encoded;
    }

    std::optional<std::vector<std::byte>> ReadContainerHeaderBytes()
    {
        auto first = ReadBytes(4);
        if (std::empty(first)) {
            return std::nullopt;  // clean EOF
        }
        if (std::size(first) < 4) {
            throw std::runtime_error("CramReader: truncated container header prefix");
        }

        std::vector<std::byte> headerBytes = first;

        // Probe explicit EOF marker when first 4 bytes match.
        if (std::ranges::equal(std::span{CRAM_EOF_MARKER}.first(4), headerBytes)) {
            const auto eofTail = ReadBytes(std::size(CRAM_EOF_MARKER) - 4);
            if (std::size(eofTail) != std::size(CRAM_EOF_MARKER) - 4) {
                throw std::runtime_error("CramReader: truncated EOF marker");
            }

            std::array<std::byte, 38> eofCandidate{};
            std::ranges::copy(headerBytes, std::begin(eofCandidate));
            std::ranges::copy(eofTail, std::begin(eofCandidate) + 4);

            // Byte 8 is an ITF-8 length nibble; only its low 4 bits are significant, so
            // mask before comparing (matches htslib cram_check_EOF and the ITF-8 decoder).
            eofCandidate[8] &= std::byte{0x0f};

            if (std::ranges::equal(CRAM_EOF_MARKER, eofCandidate)) {
                atEof = true;
                return std::nullopt;
            }

            // Not actually EOF marker; rewind and continue parsing as a container
            // header.
            file.clear();
            file.seekg(-static_cast<std::streamoff>(std::size(eofTail)), std::ios::cur);
            if (!file) {
                throw std::runtime_error("CramReader: failed to rewind after EOF probe");
            }
        }

        AppendBytes(headerBytes, ReadVarIntBytes(Itf8EncodedLength,
                                                 "container ref_seq_id"));  // ref_seq_id
        AppendBytes(headerBytes, ReadVarIntBytes(Itf8EncodedLength,
                                                 "container start_pos"));  // start_pos
        AppendBytes(headerBytes, ReadVarIntBytes(Itf8EncodedLength,
                                                 "container alignment_span"));  // alignment_span
        AppendBytes(headerBytes, ReadVarIntBytes(Itf8EncodedLength,
                                                 "container num_records"));  // num_records
        AppendBytes(headerBytes, ReadVarIntBytes(Ltf8EncodedLength,
                                                 "container record_counter"));  // record_counter
        AppendBytes(headerBytes, ReadVarIntBytes(Ltf8EncodedLength, "container bases"));  // bases
        AppendBytes(headerBytes, ReadVarIntBytes(Itf8EncodedLength,
                                                 "container num_blocks"));  // num_blocks

        const auto landmarkCountEncoded =
            ReadVarIntBytes(Itf8EncodedLength, "container landmark_count");
        AppendBytes(headerBytes, landmarkCountEncoded);
        const auto landmarkCount = DecodeEncodedItf8(landmarkCountEncoded);
        if (landmarkCount < 0) {
            throw std::runtime_error("CramReader: negative landmark count in container header");
        }
        for (std::int32_t i = 0; i < landmarkCount; ++i) {
            AppendBytes(headerBytes, ReadVarIntBytes(Itf8EncodedLength, "container landmark"));
        }

        const auto crc = ReadExactBytes(4, "container header CRC32");
        AppendBytes(headerBytes, crc);

        return headerBytes;
    }

    CramContainerHeader ParseContainerHeaderFromStream()
    {
        const auto headerBytesOpt = ReadContainerHeaderBytes();
        if (!headerBytesOpt) {
            throw std::runtime_error("CramReader: missing container header");
        }
        std::size_t bytesRead = 0;
        const auto parsedHeader = ParseContainerHeader(*headerBytesOpt, bytesRead);
        if (bytesRead != std::size(*headerBytesOpt)) {
            throw std::runtime_error("CramReader: failed to consume full container header");
        }
        return parsedHeader;
    }

    void ReadHeaderContainer()
    {
        const auto containerHeader = ParseContainerHeaderFromStream();
        if (containerHeader.Length < 0) {
            throw std::runtime_error("CramReader: negative header container length");
        }

        const auto blockData = ReadExactBytes(static_cast<std::size_t>(containerHeader.Length),
                                              "header container payload");
        std::size_t blockBytesRead = 0;
        auto block = ParseBlock(blockData, blockBytesRead);
        if (block.ContentType != CramBlockContentType::FILE_HEADER) {
            throw std::runtime_error("CramReader: first header block is not FILE_HEADER");
        }

        DecompressCramBlock(block, gzipContext.Decompressor.get());

        // The block data contains: int32 header_length + header_text
        if (std::size(block.Data) < 4) {
            throw std::runtime_error("CramReader: empty header block");
        }
        const std::int32_t textLen = ReadI32LE(block.Data.data());
        if (textLen < 0 || static_cast<std::size_t>(textLen) + 4 > std::size(block.Data)) {
            throw std::runtime_error("CramReader: invalid header text length");
        }

        const std::string headerText(reinterpret_cast<const char*>(block.Data.data() + 4), textLen);
        auto result = SamHeader::FromText(headerText);
        if (!result) {
            throw std::runtime_error{
                std::format("CramReader: failed to parse SAM header: {}", result.error())};
        }
        header = std::move(*result);
    }

    std::vector<BamRecord> DecodeContainerRecords(
        const CramContainerHeader& containerHeader, std::span<const std::byte> blockData,
        const std::optional<SliceQueryCandidate>& targetSlice = std::nullopt)
    {
        std::vector<BamRecord> decodedRecords;
        if (containerHeader.NumBlocks == 0) {
            return decodedRecords;
        }

        const auto container = ParseContainer(containerHeader, blockData);
        const auto& compressionHeader = container.CompressionHeader;
        const auto substitutionMatrix =
            BuildSubstitutionMatrixLookup(compressionHeader.PreservationMap.SubstitutionMatrix);

        std::vector<std::unique_ptr<CramCodec>> codecStorage;
        codecStorage.reserve(std::size(compressionHeader.DataSeriesEncodings));
        std::vector<std::pair<CramDataSeries, CramCodec*>> flatCodecs;
        flatCodecs.reserve(std::size(compressionHeader.DataSeriesEncodings));
        for (const auto& [key, desc] : compressionHeader.DataSeriesEncodings) {
            codecStorage.push_back(CreateCodec(desc));
            flatCodecs.emplace_back(key, codecStorage.back().get());
        }
        std::ranges::sort(flatCodecs, {}, &std::pair<CramDataSeries, CramCodec*>::first);

        const auto tagDictionary =
            ParseTagDictionary(compressionHeader.PreservationMap.TagIdsDictionary);

        std::vector<std::unique_ptr<CramCodec>> tagCodecStorage;
        tagCodecStorage.reserve(std::size(compressionHeader.TagEncodings));
        std::vector<std::pair<std::int32_t, CramCodec*>> flatTagCodecs;
        flatTagCodecs.reserve(std::size(compressionHeader.TagEncodings));
        for (const auto& [tagKey, desc] : compressionHeader.TagEncodings) {
            tagCodecStorage.push_back(CreateCodec(desc));
            flatTagCodecs.emplace_back(tagKey, tagCodecStorage.back().get());
        }
        std::ranges::sort(flatTagCodecs, {}, &std::pair<std::int32_t, CramCodec*>::first);

        std::int64_t nextComputedSliceOffset = 0;
        if (std::empty(containerHeader.Landmarks)) {
            const auto compHdrData = SerializeCompressionHeader(container.CompressionHeader);
            CramBlock compHdrBlock;
            compHdrBlock.Method = CramBlockMethod::RAW;
            compHdrBlock.ContentType = CramBlockContentType::COMPRESSION_HEADER;
            compHdrBlock.ContentId = 0;
            compHdrBlock.RawSize = static_cast<std::int32_t>(std::size(compHdrData));
            compHdrBlock.CompressedSize = compHdrBlock.RawSize;
            compHdrBlock.Data = compHdrData;
            nextComputedSliceOffset =
                static_cast<std::int64_t>(std::size(SerializeBlock(compHdrBlock)));
        }

        for (std::size_t sliceIndex = 0; sliceIndex < std::size(container.Slices); ++sliceIndex) {
            const auto& slice = container.Slices[sliceIndex];
            const auto& sliceHeader = slice.Header;
            if (sliceHeader.NumRecords < 0) {
                throw std::runtime_error("CramReader: negative slice record count");
            }

            std::int64_t sliceOffset = nextComputedSliceOffset;
            std::int64_t sliceSize = static_cast<std::int64_t>(SerializedSliceSize(slice));
            if (!std::empty(containerHeader.Landmarks)) {
                sliceOffset = containerHeader.Landmarks.at(sliceIndex);
                std::int64_t nextLandmark = static_cast<std::int64_t>(containerHeader.Length);
                if (sliceIndex + 1 < std::size(containerHeader.Landmarks)) {
                    nextLandmark =
                        static_cast<std::int64_t>(containerHeader.Landmarks.at(sliceIndex + 1));
                }
                sliceSize = nextLandmark - sliceOffset;
            }
            if (sliceSize <= 0) {
                throw std::runtime_error("CramReader: non-positive slice size");
            }
            nextComputedSliceOffset = sliceOffset + sliceSize;

            const bool shouldDecodeSlice =
                !targetSlice ||
                ((targetSlice->SliceOffset == sliceOffset) &&
                 ((targetSlice->SliceSize <= 0) || (targetSlice->SliceSize == sliceSize)));
            if (!shouldDecodeSlice) {
                continue;
            }

            std::vector<CramBlock> sliceBlocks;
            sliceBlocks.reserve(1 + std::size(slice.ExternalBlocks));
            sliceBlocks.push_back(slice.CoreBlock);
            sliceBlocks.insert(std::end(sliceBlocks), std::begin(slice.ExternalBlocks),
                               std::end(slice.ExternalBlocks));
            DecompressBlocks(sliceBlocks);

            CramBlock coreBlock;
            bool sawCoreBlock = false;
            CramExternalBlockStore extStore;
            for (auto& block : sliceBlocks) {
                if (block.ContentType == CramBlockContentType::CORE_DATA) {
                    if (sawCoreBlock) {
                        throw std::runtime_error("CramReader: multiple CORE_DATA blocks in slice");
                    }
                    coreBlock = std::move(block);
                    sawCoreBlock = true;
                } else if (block.ContentType == CramBlockContentType::EXTERNAL_DATA) {
                    extStore.AddBlock(block.ContentId, block.Data);
                } else {
                    throw std::runtime_error("CramReader: unexpected block type in slice payload");
                }
            }
            if (!sawCoreBlock) {
                throw std::runtime_error("CramReader: slice missing CORE_DATA block");
            }

            std::span<const std::byte> embeddedReference{};
            if (sliceHeader.EmbeddedRefBlockId >= 0) {
                embeddedReference = extStore.GetBlockData(sliceHeader.EmbeddedRefBlockId);
            }
            ValidateSliceReferenceMd5(sliceHeader,
                                      compressionHeader.PreservationMap.ReferenceRequired,
                                      embeddedReference, externalReferenceById, header);

            CramBitReader coreReader{coreBlock.Data};
            std::int64_t previousApos{sliceHeader.AlignmentStart};
            std::vector<BamRecord> sliceRecords;
            std::vector<std::int32_t> downstreamMateDistances;
            sliceRecords.reserve(sliceHeader.NumRecords);
            downstreamMateDistances.reserve(sliceHeader.NumRecords);

            std::int32_t downstreamMateDistance{-1};
            DecodeRecordContext decodeCtx{
                .CoreReader = coreReader,
                .ExtStore = extStore,
                .Codecs = flatCodecs,
                .TagCodecs = flatTagCodecs,
                .TagDictionary = tagDictionary,
                .SliceHeader = sliceHeader,
                .SubstitutionMatrix = substitutionMatrix,
                .EmbeddedReference = embeddedReference,
                .PreviousApos = previousApos,
                .DownstreamMateDistance = downstreamMateDistance,
                .ApDelta = compressionHeader.PreservationMap.ApDelta,
                .ReadNamesIncluded = compressionHeader.PreservationMap.ReadNamesIncluded,
                .ExternalReferenceById = externalReferenceById,
            };

            for (std::int32_t ri = 0; ri < sliceHeader.NumRecords; ++ri) {
                BamRecord record;
                downstreamMateDistance = -1;
                DecodeRecord(record, decodeCtx);
                sliceRecords.push_back(std::move(record));
                downstreamMateDistances.push_back(downstreamMateDistance);
            }

            ResolveDownstreamMates(sliceRecords, downstreamMateDistances);
            decodedRecords.insert(std::end(decodedRecords),
                                  std::make_move_iterator(std::begin(sliceRecords)),
                                  std::make_move_iterator(std::end(sliceRecords)));
        }

        return decodedRecords;
    }

    void DecompressBlocks(std::vector<CramBlock>& blocks)
    {
        if (std::empty(blocks)) {
            return;
        }
        if (std::size(blocks) >
            static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max())) {
            throw std::runtime_error("CramReader: too many blocks to decompress in one slice");
        }

        const std::size_t workers = config.DecompressionWorkers;
        constexpr std::size_t MIN_BLOCKS_FOR_PARALLEL = 8;
        if (!decompressionPool || workers <= 1 || std::size(blocks) < MIN_BLOCKS_FOR_PARALLEL) {
            for (auto& block : blocks) {
                DecompressCramBlock(block, gzipContext.Decompressor.get());
            }
            return;
        }

        std::atomic<std::int32_t> nextBlock{0};
        const std::int32_t totalBlocks = static_cast<std::int32_t>(std::size(blocks));
        const std::int32_t taskCount = std::ranges::min(workers, std::size(blocks));
        Parallel::Dispatch(
            decompressionPool,
            MakeWorkStealingTask(&nextBlock, totalBlocks,
                                 [&blocks](std::int32_t i) {
                                     auto& block = blocks[static_cast<std::size_t>(i)];
                                     DecompressCramBlock(block, ThreadLocalGzipDecompressor());
                                 }),
            taskCount);
    }

    std::vector<BamRecord> LoadIndexedSliceRecords(const SliceQueryCandidate& candidate)
    {
        if (candidate.ContainerOffset < 0 || candidate.SliceOffset < 0) {
            throw std::runtime_error("CramReader: CRAI entry has negative container/slice offset");
        }

        file.clear();
        file.seekg(static_cast<std::streamoff>(candidate.ContainerOffset), std::ios::beg);
        if (!file) {
            throw std::runtime_error("CramReader: failed to seek to CRAI container offset");
        }

        const auto containerHeader = ParseContainerHeaderFromStream();
        if (containerHeader.Length < 0) {
            throw std::runtime_error("CramReader: negative container length");
        }
        if (containerHeader.NumBlocks < 0) {
            throw std::runtime_error("CramReader: negative container block count");
        }

        const auto blockData =
            ReadExactBytes(static_cast<std::size_t>(containerHeader.Length), "container payload");
        return DecodeContainerRecords(containerHeader, blockData, candidate);
    }

    bool LoadNextContainer()
    {
        pendingRecords.clear();
        pendingIdx = 0;

        while (!atEof && std::empty(pendingRecords)) {
            const auto headerBytesOpt = ReadContainerHeaderBytes();
            if (!headerBytesOpt) {
                atEof = true;
                break;
            }

            std::size_t headerBytesRead = 0;
            const auto containerHeader = ParseContainerHeader(*headerBytesOpt, headerBytesRead);
            if (headerBytesRead != std::size(*headerBytesOpt)) {
                throw std::runtime_error("CramReader: failed to consume full container header");
            }
            if (containerHeader.Length < 0) {
                throw std::runtime_error("CramReader: negative container length");
            }
            if (containerHeader.NumBlocks < 0) {
                throw std::runtime_error("CramReader: negative container block count");
            }

            const auto blockData = ReadExactBytes(static_cast<std::size_t>(containerHeader.Length),
                                                  "container payload");
            pendingRecords = DecodeContainerRecords(containerHeader, blockData);
        }

        return !std::empty(pendingRecords);
    }

    void DecodeRecord(BamRecord& record, DecodeRecordContext& ctx)
    {
        // BF - BAM bit flags
        std::uint16_t bamFlags = 0;
        if (auto* codec = LookupCodec(ctx, CramDataSeries::BF)) {
            bamFlags = codec->DecodeInt(ctx.CoreReader, ctx.ExtStore);
        }

        // CF - CRAM bit flags
        std::int32_t cramFlags = 0;
        if (auto* codec = LookupCodec(ctx, CramDataSeries::CF)) {
            cramFlags = codec->DecodeInt(ctx.CoreReader, ctx.ExtStore);
        }

        // RI - reference sequence id (only in multi-ref slices)
        std::int32_t refId = ctx.SliceHeader.RefSeqId;
        if (ctx.SliceHeader.RefSeqId == -2) {
            if (auto* codec = LookupCodec(ctx, CramDataSeries::RI)) {
                refId = codec->DecodeInt(ctx.CoreReader, ctx.ExtStore);
            }
        }

        // RL - read length
        std::int32_t readLength = 0;
        if (auto* codec = LookupCodec(ctx, CramDataSeries::RL)) {
            readLength = codec->DecodeInt(ctx.CoreReader, ctx.ExtStore);
        }

        // AP - alignment position
        std::int32_t pos = 0;
        if (auto* codec = LookupCodec(ctx, CramDataSeries::AP)) {
            const auto apValue =
                static_cast<std::int64_t>(codec->DecodeInt(ctx.CoreReader, ctx.ExtStore));
            if (ctx.ApDelta) {
                const auto absoluteApos = ctx.PreviousApos + apValue;
                ctx.PreviousApos = absoluteApos;
                pos = absoluteApos - 1;
            } else {
                ctx.PreviousApos = apValue;
                pos = apValue - 1;
            }
        }

        // RG - read group
        std::int32_t readGroup = -1;
        if (auto* codec = LookupCodec(ctx, CramDataSeries::RG)) {
            readGroup = codec->DecodeInt(ctx.CoreReader, ctx.ExtStore);
        }

        // RN - read name (if preserved)
        std::string readName;
        if (ctx.ReadNamesIncluded) {
            if (auto* codec = LookupCodec(ctx, CramDataSeries::RN)) {
                DecodeReadName(readName, *codec, ctx.CoreReader, ctx.ExtStore);
            }
        }

        // MF, NS, NP, TS - mate info (if detached)
        std::int32_t nextRefId = -1;
        std::int32_t nextPos = -1;
        std::int32_t tlen = 0;
        ctx.DownstreamMateDistance = -1;
        if (cramFlags & CRAM_FLAG_DETACHED) {
            if (auto* codec = LookupCodec(ctx, CramDataSeries::MF)) {
                codec->DecodeInt(ctx.CoreReader, ctx.ExtStore);  // mate flags
            }
            if (auto* codec = LookupCodec(ctx, CramDataSeries::NS)) {
                nextRefId = codec->DecodeInt(ctx.CoreReader, ctx.ExtStore);
            }
            if (auto* codec = LookupCodec(ctx, CramDataSeries::NP)) {
                nextPos = codec->DecodeInt(ctx.CoreReader, ctx.ExtStore) - 1;  // 0-based
            }
            if (auto* codec = LookupCodec(ctx, CramDataSeries::TS)) {
                tlen = codec->DecodeInt(ctx.CoreReader, ctx.ExtStore);
            }

            // RN after mate info if not preserved globally
            if (!ctx.ReadNamesIncluded) {
                if (auto* codec = LookupCodec(ctx, CramDataSeries::RN)) {
                    DecodeReadName(readName, *codec, ctx.CoreReader, ctx.ExtStore);
                }
            }
        } else if (cramFlags & CRAM_FLAG_HAS_MATE_DOWNSTREAM) {
            if (auto* codec = LookupCodec(ctx, CramDataSeries::NF)) {
                ctx.DownstreamMateDistance = codec->DecodeInt(ctx.CoreReader, ctx.ExtStore);
            }
        }

        // TL - tag IDs dictionary index
        std::int32_t tagListIndex{0};
        if (auto* codec = LookupCodec(ctx, CramDataSeries::TL)) {
            tagListIndex = codec->DecodeInt(ctx.CoreReader, ctx.ExtStore);
        }

        // Decode sequence data
        std::string sequence;
        std::vector<std::uint8_t> qualities;
        std::vector<CigarOp> cigar;
        std::uint8_t mapQ = 0;

        if ((bamFlags & 0x4) == 0) {
            // Mapped read
            // FN - number of read features
            std::int32_t numFeatures = 0;
            if (auto* codec = LookupCodec(ctx, CramDataSeries::FN)) {
                numFeatures = codec->DecodeInt(ctx.CoreReader, ctx.ExtStore);
            }

            std::vector<DecodedFeature> features;
            features.reserve(numFeatures);

            // Read features
            std::int32_t prevFeaturePos = 0;
            for (std::int32_t fi = 0; fi < numFeatures; ++fi) {
                DecodedFeature feature;
                if (auto* codec = LookupCodec(ctx, CramDataSeries::FC)) {
                    feature.Code = static_cast<char>(
                        static_cast<std::uint8_t>(codec->DecodeByte(ctx.CoreReader, ctx.ExtStore)));
                }

                if (auto* codec = LookupCodec(ctx, CramDataSeries::FP)) {
                    const auto delta = codec->DecodeInt(ctx.CoreReader, ctx.ExtStore);
                    prevFeaturePos += delta;
                    feature.ReadPos = prevFeaturePos;
                }

                switch (feature.Code) {
                    case 'B':
                        // Read base: one base (BA) plus one quality (QS); consumes an M.
                        if (auto* baseCodec = LookupCodec(ctx, CramDataSeries::BA)) {
                            feature.Data.push_back(
                                baseCodec->DecodeByte(ctx.CoreReader, ctx.ExtStore));
                            feature.IsSubstitutionCode = false;
                        }
                        if (auto* qsCodec = LookupCodec(ctx, CramDataSeries::QS)) {
                            feature.Quality = static_cast<std::uint8_t>(
                                qsCodec->DecodeByte(ctx.CoreReader, ctx.ExtStore));
                        }
                        feature.Length = 1;
                        break;
                    case 'X':
                        if (auto* codec = LookupCodec(ctx, CramDataSeries::BS)) {
                            feature.Data.push_back(codec->DecodeByte(ctx.CoreReader, ctx.ExtStore));
                            feature.IsSubstitutionCode = true;
                        }
                        feature.Length = 1;
                        break;
                    case 'I':
                        if (auto* codec = LookupCodec(ctx, CramDataSeries::IN)) {
                            feature.Data = codec->DecodeByteArray(ctx.CoreReader, ctx.ExtStore);
                            feature.Length = std::size(feature.Data);
                        }
                        break;
                    case 'D':
                        if (auto* codec = LookupCodec(ctx, CramDataSeries::DL)) {
                            feature.Length = codec->DecodeInt(ctx.CoreReader, ctx.ExtStore);
                        }
                        break;
                    case 'i':
                        if (auto* codec = LookupCodec(ctx, CramDataSeries::BA)) {
                            feature.Data.push_back(codec->DecodeByte(ctx.CoreReader, ctx.ExtStore));
                            feature.Length = 1;
                        }
                        break;
                    case 'b':
                        if (auto* codec = LookupCodec(ctx, CramDataSeries::BB)) {
                            feature.Data = codec->DecodeByteArray(ctx.CoreReader, ctx.ExtStore);
                            feature.Length = std::size(feature.Data);
                        }
                        break;
                    case 'q':
                        // Several quality values: array from QQ.
                        if (auto* codec = LookupCodec(ctx, CramDataSeries::QQ)) {
                            feature.Data = codec->DecodeByteArray(ctx.CoreReader, ctx.ExtStore);
                            feature.Length = std::size(feature.Data);
                        }
                        break;
                    case 'Q':
                        // Quality score: single byte from QS.
                        if (auto* codec = LookupCodec(ctx, CramDataSeries::QS)) {
                            feature.Data.push_back(codec->DecodeByte(ctx.CoreReader, ctx.ExtStore));
                            feature.Length = 1;
                        }
                        break;
                    case 'S':
                        if (auto* codec = LookupCodec(ctx, CramDataSeries::SC)) {
                            feature.Data = codec->DecodeByteArray(ctx.CoreReader, ctx.ExtStore);
                            feature.Length = std::size(feature.Data);
                        }
                        break;
                    case 'N':
                        if (auto* codec = LookupCodec(ctx, CramDataSeries::RS)) {
                            feature.Length = codec->DecodeInt(ctx.CoreReader, ctx.ExtStore);
                        }
                        break;
                    case 'P':
                        if (auto* codec = LookupCodec(ctx, CramDataSeries::PD)) {
                            feature.Length = codec->DecodeInt(ctx.CoreReader, ctx.ExtStore);
                        }
                        break;
                    case 'H':
                        if (auto* codec = LookupCodec(ctx, CramDataSeries::HC)) {
                            feature.Length = codec->DecodeInt(ctx.CoreReader, ctx.ExtStore);
                        }
                        break;
                    default:
                        break;
                }

                features.push_back(std::move(feature));
            }

            // MQ - mapping quality
            if (auto* codec = LookupCodec(ctx, CramDataSeries::MQ)) {
                mapQ = codec->DecodeInt(ctx.CoreReader, ctx.ExtStore);
            }

            // QS - quality scores (if stored as array)
            if (cramFlags & CRAM_FLAG_QUALITY_AS_ARRAY) {
                if (auto* qsCodec = LookupCodec(ctx, CramDataSeries::QS)) {
                    DecodeQualityArray(qualities, *qsCodec, ctx.CoreReader, ctx.ExtStore,
                                       readLength);
                }
            } else if (readLength > 0) {
                qualities.assign(readLength, 0xFF);
                for (const auto& feature : features) {
                    if (feature.ReadPos <= 0) {
                        continue;
                    }
                    const auto qualPos = static_cast<std::size_t>(feature.ReadPos - 1);
                    if (feature.Code == 'q') {
                        // Several quality values (QQ array): a run starting at ReadPos-1.
                        if (!std::empty(feature.Data)) {
                            CopyBytesToVector(feature.Data, qualities, qualPos);
                        }
                    } else if (feature.Code == 'Q') {
                        // Quality score (single QS byte) placed at ReadPos-1.
                        if (!std::empty(feature.Data) && qualPos < std::size(qualities)) {
                            qualities[qualPos] = static_cast<std::uint8_t>(feature.Data.front());
                        }
                    } else if (feature.Code == 'B') {
                        // Read base: QS byte captured alongside the base.
                        if (feature.Quality.has_value() && qualPos < std::size(qualities)) {
                            qualities[qualPos] = *feature.Quality;
                        }
                    }
                }
            }

            // Reconstruct sequence and CIGAR.
            if (readLength > 0) {
                sequence.assign(readLength, 'N');
            }

            std::int32_t currentReadPos{1};
            std::int32_t currentRefPos{pos + 1};  // 1-based reference coordinate
            for (const auto& feature : features) {
                if (feature.ReadPos > currentReadPos) {
                    const auto matchLen = feature.ReadPos - currentReadPos;
                    AppendCigarOp(cigar, CigarOpType::M, matchLen);
                    FillReferenceMatches(sequence, ctx, refId, currentReadPos, currentRefPos,
                                         matchLen);
                    currentReadPos += matchLen;
                    currentRefPos += matchLen;
                }

                switch (feature.Code) {
                    case 'b':
                        if (feature.Length > 0) {
                            WriteFeatureData(sequence, currentReadPos, feature.Data);
                            AppendCigarOp(cigar, CigarOpType::M, feature.Length);
                            currentReadPos += feature.Length;
                            currentRefPos += feature.Length;
                        }
                        break;
                    case 'B':
                    case 'X': {
                        if (currentReadPos > 0 &&
                            static_cast<std::size_t>(currentReadPos - 1) < std::size(sequence)) {
                            char base = 'N';
                            if (!std::empty(feature.Data)) {
                                if (feature.IsSubstitutionCode) {
                                    base = SubstitutionBase(
                                        ctx, refId, static_cast<std::uint8_t>(feature.Data.front()),
                                        currentRefPos);
                                } else {
                                    base = static_cast<char>(
                                        static_cast<std::uint8_t>(feature.Data.front()));
                                }
                            }
                            sequence[static_cast<std::size_t>(currentReadPos - 1)] = base;
                        }
                        AppendCigarOp(cigar, CigarOpType::M, 1);
                        currentReadPos += 1;
                        currentRefPos += 1;
                        break;
                    }
                    case 'I': {
                        const std::int32_t insLen{feature.Length};
                        if (insLen > 0) {
                            WriteFeatureData(sequence, currentReadPos, feature.Data);
                            AppendCigarOp(cigar, CigarOpType::I, insLen);
                            currentReadPos += insLen;
                        }
                        break;
                    }
                    case 'i':
                        if (!std::empty(feature.Data)) {
                            WriteFeatureData(sequence, currentReadPos, {feature.Data.data(), 1});
                        }
                        AppendCigarOp(cigar, CigarOpType::I, 1);
                        currentReadPos += 1;
                        break;
                    case 'S':
                        if (feature.Length > 0) {
                            WriteFeatureData(sequence, currentReadPos, feature.Data);
                            AppendCigarOp(cigar, CigarOpType::S, feature.Length);
                            currentReadPos += feature.Length;
                        }
                        break;
                    case 'D':
                        if (feature.Length > 0) {
                            AppendCigarOp(cigar, CigarOpType::D, feature.Length);
                            currentRefPos += feature.Length;
                        }
                        break;
                    case 'N':
                        if (feature.Length > 0) {
                            AppendCigarOp(cigar, CigarOpType::N, feature.Length);
                            currentRefPos += feature.Length;
                        }
                        break;
                    case 'H': {
                        AppendCigarOp(cigar, CigarOpType::H, FeatureLengthOrOne(feature));
                        break;
                    }
                    case 'P': {
                        AppendCigarOp(cigar, CigarOpType::P, FeatureLengthOrOne(feature));
                        break;
                    }
                    default:
                        break;
                }
            }

            if (readLength >= currentReadPos) {
                const auto tailLen = readLength - currentReadPos + 1;
                AppendCigarOp(cigar, CigarOpType::M, tailLen);
                FillReferenceMatches(sequence, ctx, refId, currentReadPos, currentRefPos, tailLen);
            }
            if (std::empty(cigar) && readLength > 0) {
                AppendCigarOp(cigar, CigarOpType::M, readLength);
            }
        } else {
            // Unmapped read
            sequence.resize(readLength);
            if (auto* baCodec = LookupCodec(ctx, CramDataSeries::BA)) {
                for (std::int32_t bi = 0; bi < readLength; ++bi) {
                    sequence[bi] = static_cast<char>(static_cast<std::uint8_t>(
                        baCodec->DecodeByte(ctx.CoreReader, ctx.ExtStore)));
                }
            }

            // QS quality scores
            if (cramFlags & CRAM_FLAG_QUALITY_AS_ARRAY) {
                if (auto* qsCodec = LookupCodec(ctx, CramDataSeries::QS)) {
                    DecodeQualityArray(qualities, *qsCodec, ctx.CoreReader, ctx.ExtStore,
                                       readLength);
                }
            }
        }

        if (cramFlags & CRAM_FLAG_SEQUENCE_OMITTED) {
            sequence.clear();
            qualities.clear();
        }

        // Build the BamRecord
        record.Name(std::move(readName))
            .Flag(bamFlags)
            .RefId(refId)
            .Pos(pos)
            .MapQ(mapQ)
            .Cigar(std::move(cigar))
            .NextRefId(nextRefId)
            .NextPos(nextPos)
            .Tlen(tlen)
            .Sequence(std::move(sequence))
            .Qualities(std::move(qualities));

        if (tagListIndex >= 0 &&
            tagListIndex < static_cast<std::int32_t>(std::size(ctx.TagDictionary))) {
            for (const auto& tag : ctx.TagDictionary[tagListIndex]) {
                const auto contentId = TagContentId(tag.Tag1, tag.Tag2, tag.Type);
                if (auto* codec = LookupTagCodec(ctx, contentId)) {
                    const auto payload =
                        DecodeTagPayload(tag.Type, *codec, ctx.CoreReader, ctx.ExtStore);
                    record.MutableTags().Set(TagKey{tag.Tag1, tag.Tag2},
                                             DecodeTagValueFromBamPayload(tag.Type, payload));
                }
            }
        }

        if (readGroup >= 0 &&
            readGroup < static_cast<std::int32_t>(std::size(header.ReadGroups()))) {
            record.MutableTags().Set(RG_TAG, std::string{header.ReadGroups()[readGroup].Id()});
        }
    }
};

// ---------------------------------------------------------------------------
// CramReader
// ---------------------------------------------------------------------------

namespace {

std::optional<BamRecord> TakePendingRecord(std::vector<BamRecord>& pendingRecords,
                                           std::size_t& pendingIdx)
{
    if (pendingIdx >= std::size(pendingRecords)) {
        return std::nullopt;
    }

    return std::move(pendingRecords[pendingIdx++]);
}

}  // namespace

CramReader::CramReader(const std::filesystem::path& path, const CramReaderConfig& config)
    : impl_{std::make_unique<Impl>()}
{
    impl_->path = path;
    impl_->config = config;
    impl_->Open();
}

CramReader::~CramReader() = default;
CramReader::CramReader(CramReader&&) noexcept = default;
CramReader& CramReader::operator=(CramReader&&) noexcept = default;

const SamHeader& CramReader::Header() const { return impl_->header; }

std::optional<BamRecord> CramReader::ReadRecord()
{
    if (auto record{TakePendingRecord(impl_->pendingRecords, impl_->pendingIdx)}) {
        return record;
    }

    if (impl_->atEof || !impl_->LoadNextContainer()) {
        return std::nullopt;
    }

    return TakePendingRecord(impl_->pendingRecords, impl_->pendingIdx);
}

std::optional<RawRecord> CramReader::ReadRawRecord()
{
    auto record = ReadRecord();
    if (!record) {
        return std::nullopt;
    }
    return ToRawRecord(*record);
}

std::vector<BamRecord> CramReader::Query(const CraiIndex& index, std::int32_t refId,
                                         std::int32_t beg, std::int32_t end)
{
    if ((refId < -1) || (refId >= 0 && end <= beg)) {
        return {};
    }

    const std::vector<SliceQueryCandidate> candidates =
        BuildSliceCandidates(index, refId, beg, end);
    if (std::empty(candidates)) {
        return {};
    }

    // Query execution uses an isolated reader instance so sequential iteration
    // state on this reader remains unchanged.
    CramReader queryReader{impl_->path, impl_->config};
    std::vector<BamRecord> queriedRecords;

    for (const SliceQueryCandidate& candidate : candidates) {
        std::vector<BamRecord> decoded = queryReader.impl_->LoadIndexedSliceRecords(candidate);
        for (BamRecord& record : decoded) {
            if (KeepRecordForQuery(record, refId, beg, end)) {
                queriedRecords.push_back(std::move(record));
            }
        }
    }

    return queriedRecords;
}

std::vector<RawRecord> CramReader::QueryRaw(const CraiIndex& index, std::int32_t refId,
                                            std::int32_t beg, std::int32_t end)
{
    const std::vector<BamRecord> records = Query(index, refId, beg, end);
    std::vector<RawRecord> rawRecords;
    rawRecords.reserve(std::size(records));
    for (const BamRecord& record : records) {
        rawRecords.push_back(ToRawRecord(record));
    }
    return rawRecords;
}

// ---------------------------------------------------------------------------
// RecordRange
// ---------------------------------------------------------------------------

CramReader::RecordRange::RecordRange(CramReader* reader) : reader_{reader} {}

CramReader::RecordRange::Iterator CramReader::RecordRange::begin() { return Iterator{reader_}; }

CramReader::RecordRange::Iterator CramReader::RecordRange::end() { return Iterator{}; }

CramReader::RecordRange CramReader::Records() { return RecordRange{this}; }

// --- Iterator ---

CramReader::RecordRange::Iterator::Iterator() = default;

CramReader::RecordRange::Iterator::Iterator(CramReader* reader) : reader_{reader}
{
    detail::AdvanceReaderIterator(reader_, current_);
}

const BamRecord& CramReader::RecordRange::Iterator::operator*() const { return *current_; }

const BamRecord* CramReader::RecordRange::Iterator::operator->() const { return &*current_; }

CramReader::RecordRange::Iterator& CramReader::RecordRange::Iterator::operator++()
{
    detail::AdvanceReaderIterator(reader_, current_);
    return *this;
}

void CramReader::RecordRange::Iterator::operator++(int) { ++(*this); }

bool CramReader::RecordRange::Iterator::operator==(const Iterator& other) const
{
    return reader_ == other.reader_;
}

// ---------------------------------------------------------------------------
// RawRecordRange
// ---------------------------------------------------------------------------

CramReader::RawRecordRange::RawRecordRange(CramReader* reader) : reader_{reader} {}

CramReader::RawRecordRange::Iterator CramReader::RawRecordRange::begin()
{
    return Iterator{reader_};
}

CramReader::RawRecordRange::Iterator CramReader::RawRecordRange::end() { return Iterator{}; }

CramReader::RawRecordRange CramReader::RawRecords() { return RawRecordRange{this}; }

CramReader::RawRecordRange::Iterator::Iterator() = default;

CramReader::RawRecordRange::Iterator::Iterator(CramReader* reader) : reader_{reader}
{
    detail::AdvanceReaderIterator(reader_, current_, &CramReader::ReadRawRecord);
}

const RawRecord& CramReader::RawRecordRange::Iterator::operator*() const { return *current_; }

const RawRecord* CramReader::RawRecordRange::Iterator::operator->() const { return &*current_; }

CramReader::RawRecordRange::Iterator& CramReader::RawRecordRange::Iterator::operator++()
{
    detail::AdvanceReaderIterator(reader_, current_, &CramReader::ReadRawRecord);
    return *this;
}

void CramReader::RawRecordRange::Iterator::operator++(int) { ++(*this); }

bool CramReader::RawRecordRange::Iterator::operator==(const Iterator& other) const
{
    return reader_ == other.reader_;
}

}  // namespace Samoa
}  // namespace PacBio
