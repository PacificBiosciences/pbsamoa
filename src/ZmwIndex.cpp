#include <pbsamoa/index/ZmwIndex.hpp>

#include "BinaryUtils.hpp"
#include "PathUtils.hpp"
#include "ZmiInternal.hpp"

#include <pbsamoa/core/Bgzf.hpp>

#include <array>
#include <bit>
#include <format>
#include <limits>
#include <optional>
#include <ranges>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <cstddef>
#include <cstdint>
#include <cstring>

namespace PacBio {
namespace Samoa {

ZmwIndex::ZmwIndex(ZmwIndex&& other) noexcept
    : rgIds_{std::move(other.rgIds_)}
    , zmws_{std::move(other.zmws_)}
    , offsets_{std::move(other.offsets_)}
    , queryLengths_{std::move(other.queryLengths_)}
    , hasBarcodeData_{other.hasBarcodeData_}
    , barcodeOffsets_{std::move(other.barcodeOffsets_)}
    , zmwIndex_{std::move(other.zmwIndex_)}
    , identityIndex_{std::move(other.identityIndex_)}
    , cachedNumZmws_{other.cachedNumZmws_}
{
    // indexOnce_ is left default-constructed (unset).
    // If the source had already built its index, the moved maps are valid
    // and BuildIndex will simply re-detect that via the non-empty maps.
}

ZmwIndex& ZmwIndex::operator=(ZmwIndex&& other) noexcept
{
    if (this != &other) {
        rgIds_ = std::move(other.rgIds_);
        zmws_ = std::move(other.zmws_);
        offsets_ = std::move(other.offsets_);
        queryLengths_ = std::move(other.queryLengths_);
        hasBarcodeData_ = other.hasBarcodeData_;
        barcodeOffsets_ = std::move(other.barcodeOffsets_);
        zmwIndex_ = std::move(other.zmwIndex_);
        identityIndex_ = std::move(other.identityIndex_);
        cachedNumZmws_ = other.cachedNumZmws_;
        // indexOnce_ cannot be moved; leave as default (unset).
        // call_once will rebuild if needed, but moved maps are already populated.
    }
    return *this;
}

namespace {

// PBI format constants
constexpr std::size_t PBI_HEADER_SIZE{32};
constexpr std::array<char, 4> PBI_MAGIC{'P', 'B', 'I', '\1'};

std::size_t CheckedMul(std::size_t lhs, std::size_t rhs, std::string_view what)
{
    if ((lhs != 0) && (rhs > (std::numeric_limits<std::size_t>::max() / lhs))) {
        throw std::runtime_error{std::format("{} size overflows", what)};
    }
    return lhs * rhs;
}

class BgzfDataReader
{
public:
    explicit BgzfDataReader(const std::filesystem::path& path) : reader_{path} {}

    void ReadExactly(std::span<std::byte> output, std::string_view error)
    {
        while (!std::empty(output)) {
            if (!Refill()) {
                throw std::runtime_error{std::string{error}};
            }
            const std::size_t count{std::min(std::size(output), size_ - pos_)};
            std::memcpy(std::data(output), std::data(block_) + pos_, count);
            output = output.subspan(count);
            pos_ += count;
        }
    }

    void SkipExactly(std::size_t count, std::string_view error)
    {
        while (count != 0) {
            if (!Refill()) {
                throw std::runtime_error{std::string{error}};
            }
            const std::size_t step{std::min(count, size_ - pos_)};
            pos_ += step;
            count -= step;
        }
    }

private:
    bool Refill()
    {
        if (pos_ != size_) {
            return true;
        }
        const std::optional<std::size_t> bytesRead{reader_.ReadBlock(block_)};
        if (!bytesRead) {
            throw std::runtime_error{"Failed to decompress BGZF index data"};
        }
        pos_ = 0;
        size_ = *bytesRead;
        return size_ != 0;
    }

    BgzfReader reader_;
    std::array<std::byte, MAX_DECOMPRESSED_BLOCK_SIZE> block_{};
    std::size_t pos_{0};
    std::size_t size_{0};
};

std::vector<std::byte> ReadAllBgzfData(BgzfReader& reader)
{
    std::vector<std::byte> data;
    std::array<std::byte, MAX_DECOMPRESSED_BLOCK_SIZE> blockBuf{};

    while (true) {
        const std::optional<std::size_t> bytesRead{reader.ReadBlock(blockBuf)};
        if (!bytesRead || (*bytesRead == 0)) {
            break;
        }
        data.insert(std::ranges::end(data), std::ranges::begin(blockBuf),
                    std::ranges::begin(blockBuf) + static_cast<std::ptrdiff_t>(*bytesRead));
    }

    return data;
}

template <typename T>
void ReadLittleEndianColumn(std::vector<T>& output, BgzfDataReader& reader, std::uint32_t count,
                            std::string_view error)
{
    static_assert(std::is_integral_v<T>);
    output.resize(count);
    reader.ReadExactly(std::as_writable_bytes(std::span{output}), error);
    if constexpr (std::endian::native == std::endian::big) {
        for (T& value : output) {
            using UnsignedT = std::make_unsigned_t<T>;
            const UnsignedT bits{std::bit_cast<UnsignedT>(value)};
            value = std::bit_cast<T>(std::byteswap(bits));
        }
    }
}

bool IdentityLess(const ZmwIdentity lhs, const ZmwIdentity rhs)
{
    return (lhs.rgId < rhs.rgId) || ((lhs.rgId == rhs.rgId) && (lhs.zmw < rhs.zmw));
}

std::uint64_t PartitionBoundary(std::uint64_t count, std::uint32_t boundary,
                                std::uint32_t partitions)
{
    const std::uint64_t quotient{count / partitions};
    const std::uint64_t remainder{count % partitions};
    return (quotient * boundary) + ((remainder * boundary) / partitions);
}

template <typename IndexMap, typename Key>
const typename IndexMap::mapped_type* LookupIndices(const IndexMap& index, const Key& key)
{
    if (const auto it{index.find(key)}; it != index.end()) {
        return &it->second;
    }
    return nullptr;
}

void AppendOffsets(std::vector<std::int64_t>& output, std::span<const std::ptrdiff_t> indices,
                   std::span<const std::int64_t> offsets)
{
    output.reserve(std::size(output) + std::size(indices));
    for (const std::ptrdiff_t idx : indices) {
        output.push_back(offsets[idx]);
    }
}

template <typename IndexMap, typename Key>
void AppendLookupOffsets(std::vector<std::int64_t>& output, const IndexMap& index, const Key& key,
                         std::span<const std::int64_t> offsets)
{
    if (const auto* indices{LookupIndices(index, key)}; indices) {
        AppendOffsets(output, *indices, offsets);
    }
}

template <typename IndexMap, typename Key, typename Error>
std::int64_t FirstOffsetOrThrow(const IndexMap& index, const Key& key,
                                std::span<const std::int64_t> offsets, Error&& error)
{
    const auto* indices{LookupIndices(index, key)};
    if (!indices) {
        throw std::runtime_error{error()};
    }
    return offsets[indices->front()];
}

struct ZmwNotFoundError
{
    std::int32_t zmw;

    std::string operator()() const { return std::format("ZMW not found: {}", zmw); }
};

struct ZmwIdentityNotFoundError
{
    ZmwIdentity id;

    std::string operator()() const
    {
        return std::format("ZMW identity not found: rgId={} zmw={}", id.rgId, id.zmw);
    }
};

}  // namespace

ZmwIndex ZmwIndex::FromZmi(const std::filesystem::path& path)
{
    BgzfReader reader{path};
    std::vector<std::byte> data{ReadAllBgzfData(reader)};

    // Validate minimum size for header
    if (std::size(data) < detail::ZMI_HEADER_SIZE) {
        throw std::runtime_error{"ZMI file too small: missing header"};
    }

    // Validate magic "ZMI\1"
    if (std::memcmp(std::data(data), std::data(detail::ZMI_MAGIC), 4) != 0) {
        throw std::runtime_error{"ZMI file has invalid magic bytes"};
    }

    // Read entrySize from offset 8 (2 bytes LE)
    const std::uint16_t entrySize{ReadU16LE(std::data(data) + 8)};
    if (entrySize < detail::ZMI_BASE_ENTRY_SIZE) {
        throw std::runtime_error{std::format("ZMI entrySize too small: {}", entrySize)};
    }
    const std::uint32_t version{ReadU32LE(std::data(data) + 4)};
    const bool hasBarcodes{version >= detail::ZMI_VERSION};
    if (hasBarcodes && (entrySize < detail::ZMI_ENTRY_SIZE)) {
        throw std::runtime_error{"ZMI entry too small for barcode fields"};
    }
    const std::uint64_t numRecordsHeader{ReadU64LE(std::data(data) + 12)};

    // Iterate over entries starting at offset 64 with stride = entrySize
    const std::size_t dataSize{std::size(data)};
    const std::size_t payloadSize{dataSize - detail::ZMI_HEADER_SIZE};
    if ((payloadSize % entrySize) != 0) {
        throw std::runtime_error{"ZMI file has truncated entry data"};
    }
    const std::size_t numEntries{payloadSize / entrySize};
    if ((numRecordsHeader != 0U) && (numRecordsHeader != static_cast<std::uint64_t>(numEntries))) {
        throw std::runtime_error{std::format("ZMI header numRecords mismatch: header={} entries={}",
                                             numRecordsHeader, numEntries)};
    }

    ZmwIndex index;
    index.hasBarcodeData_ = hasBarcodes;
    index.rgIds_.reserve(numEntries);
    index.zmws_.reserve(numEntries);
    index.offsets_.reserve(numEntries);

    for (std::size_t i{0}; i < numEntries; ++i) {
        const std::byte* entry{std::data(data) + detail::ZMI_HEADER_SIZE + (i * entrySize)};

        const std::int32_t rgId{ReadI32LE(entry)};
        const std::int32_t zmw{ReadI32LE(entry + 4)};
        const std::int64_t virtualOffset{ReadI64LE(entry + 8)};

        index.rgIds_.push_back(rgId);
        index.zmws_.push_back(zmw);
        index.offsets_.push_back(virtualOffset);
        if (hasBarcodes) {
            const BarcodePair pair{ReadI32LE(entry + 16), ReadI32LE(entry + 20)};
            if ((pair.Forward < -1) || (pair.Forward > 65535) || (pair.Reverse < -1) ||
                (pair.Reverse > 65535)) {
                throw std::runtime_error{"ZMI barcode IDs must be in [-1, 65535]"};
            }
            index.barcodeOffsets_[pair].push_back(virtualOffset);
        }
    }

    return index;
}

ZmwIndex ZmwIndex::FromPbi(const std::filesystem::path& path)
{
    BgzfDataReader reader{path};
    std::array<std::byte, PBI_HEADER_SIZE> header{};
    reader.ReadExactly(header, "PBI file too small: missing header");

    // Validate magic "PBI\1"
    if (std::memcmp(std::data(header), std::data(PBI_MAGIC), 4) != 0) {
        throw std::runtime_error{"PBI file has invalid magic bytes"};
    }
    // Read numReads from offset 10 (uint32 LE)
    const std::uint32_t numReads{ReadU32LE(std::data(header) + 10)};

    // BasicData columns start at offset 32 (end of header).
    // Columns are stored sequentially (not interleaved):
    //   rgId:       numReads × int32
    //   qStart:     numReads × int32
    //   qEnd:       numReads × int32
    //   holeNumber: numReads × int32
    //   readQual:   numReads × float32 (skip)
    //   ctxtFlag:   numReads × uint8   (skip)
    //   fileOffset: numReads × int64

    const std::size_t readsAsSize{numReads};
    const std::size_t float32ColBytes{CheckedMul(readsAsSize, sizeof(float), "PBI float32 column")};
    const std::size_t uint8ColBytes{
        CheckedMul(readsAsSize, sizeof(std::uint8_t), "PBI uint8 column")};
    ZmwIndex index;
    constexpr std::string_view INCOMPLETE_BASIC_DATA{
        "PBI file too small: incomplete BasicData section"};

    // Column 1: rgId (numReads × int32)
    ReadLittleEndianColumn(index.rgIds_, reader, numReads, INCOMPLETE_BASIC_DATA);

    // Columns 2 and 3: qStart and qEnd, kept as the query length qEnd - qStart
    std::vector<std::int32_t> queryStarts;
    ReadLittleEndianColumn(queryStarts, reader, numReads, INCOMPLETE_BASIC_DATA);
    ReadLittleEndianColumn(index.queryLengths_, reader, numReads, INCOMPLETE_BASIC_DATA);
    for (const auto [length, start] : std::views::zip(index.queryLengths_, queryStarts)) {
        length -= start;
    }

    // Column 4: holeNumber (numReads × int32)
    ReadLittleEndianColumn(index.zmws_, reader, numReads, INCOMPLETE_BASIC_DATA);

    // Column 5: readQual (skip)
    reader.SkipExactly(float32ColBytes, INCOMPLETE_BASIC_DATA);

    // Column 6: ctxtFlag (skip)
    reader.SkipExactly(uint8ColBytes, INCOMPLETE_BASIC_DATA);

    // Column 7: fileOffset (numReads × int64)
    ReadLittleEndianColumn(index.offsets_, reader, numReads, INCOMPLETE_BASIC_DATA);

    return index;
}

ZmwIndex ZmwIndex::Open(const std::filesystem::path& bamPath)
{
    const std::filesystem::path zmiPath{SidecarPath(bamPath, ".zmi")};
    if (std::filesystem::exists(zmiPath)) {
        return FromZmi(zmiPath);
    }
    const std::filesystem::path pbiPath{SidecarPath(bamPath, ".pbi")};
    if (std::filesystem::exists(pbiPath)) {
        return FromPbi(pbiPath);
    }
    throw std::runtime_error{std::format("No ZMW index found for: {}", bamPath.string())};
}

std::uint64_t ZmwIndex::IdentityKey(std::int32_t rgId, std::int32_t zmw)
{
    return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(rgId)) << 32) |
           static_cast<std::uint64_t>(static_cast<std::uint32_t>(zmw));
}

std::uint64_t ZmwIndex::IdentityKey(ZmwIdentity id) { return IdentityKey(id.rgId, id.zmw); }

void ZmwIndex::BuildIndex() const { std::call_once(indexOnce_, &ZmwIndex::BuildIndexImpl, this); }

void ZmwIndex::BuildIndexImpl() const
{
    // Skip if already populated (e.g., moved from a built index).
    if (!std::empty(zmwIndex_)) {
        return;
    }
    const std::ptrdiff_t n{std::ssize(zmws_)};
    zmwIndex_.reserve(std::size(zmws_));
    identityIndex_.reserve(std::size(zmws_));
    for (std::ptrdiff_t i{0}; i < n; ++i) {
        zmwIndex_[zmws_[i]].push_back(i);
        const std::uint64_t key{IdentityKey(rgIds_[i], zmws_[i])};
        identityIndex_[key].push_back(i);
    }
    cachedNumZmws_ = std::size(identityIndex_);
}

std::vector<std::int64_t> ZmwIndex::Find(std::int32_t zmw) const
{
    BuildIndex();

    std::vector<std::int64_t> result;
    AppendLookupOffsets(result, zmwIndex_, zmw, offsets_);
    return result;
}

bool ZmwIndex::HasBarcodeData() const { return hasBarcodeData_; }

std::span<const std::int64_t> ZmwIndex::FindBarcodes(const BarcodePair barcodes) const
{
    if (!hasBarcodeData_) {
        throw std::runtime_error{"Index has no barcode data; rebuild the ZMI index"};
    }
    const auto found{barcodeOffsets_.find(barcodes)};
    return (found == std::end(barcodeOffsets_)) ? std::span<const std::int64_t>{}
                                                : std::span<const std::int64_t>{found->second};
}

std::vector<BarcodePair> ZmwIndex::BarcodePairs() const
{
    if (!hasBarcodeData_) {
        throw std::runtime_error{"Index has no barcode data; rebuild the ZMI index"};
    }
    return barcodeOffsets_ | std::views::keys | std::ranges::to<std::vector>();
}

std::vector<std::int64_t> ZmwIndex::Find(ZmwIdentity id) const
{
    BuildIndex();

    std::vector<std::int64_t> result;
    AppendLookupOffsets(result, identityIndex_, IdentityKey(id), offsets_);
    return result;
}

std::vector<std::int64_t> ZmwIndex::Find(std::span<const ZmwIdentity> ids) const
{
    BuildIndex();

    std::vector<std::int64_t> result;
    result.reserve(std::size(ids));
    for (const ZmwIdentity& id : ids) {
        AppendLookupOffsets(result, identityIndex_, IdentityKey(id), offsets_);
    }
    return result;
}

std::vector<ZmwIdentity> ZmwIndex::UniqueZmws() const
{
    std::vector<ZmwIdentity> result;
    result.reserve(std::size(zmws_));

    std::unordered_set<std::uint64_t> seen;
    seen.reserve(std::size(zmws_));
    for (std::ptrdiff_t i{0}; i < std::ssize(zmws_); ++i) {
        const ZmwIdentity id{rgIds_[i], zmws_[i]};
        const std::uint64_t key{IdentityKey(id)};
        if (seen.insert(key).second) {
            result.push_back(id);
        }
    }

    return result;
}

ZmwChunk ZmwIndex::ContiguousChunk(const std::uint32_t chunkNum,
                                   const std::uint32_t totalChunks) const
{
    if ((totalChunks == 0) || (chunkNum == 0) || (chunkNum > totalChunks)) {
        throw std::invalid_argument{
            std::format("Chunk number must be in [1, {}], got {}", totalChunks, chunkNum)};
    }
    if (std::empty(zmws_)) {
        return {};
    }

    // The garden-path BAM is ordered by (rgId, ZMW), so its unique identities are
    // simply adjacent runs. Count those with no allocation or hashing. Strictly
    // increasing run keys also prove that an identity cannot recur later.
    std::uint64_t numZmws{1};
    bool ordered{true};
    ZmwIdentity previous{rgIds_.front(), zmws_.front()};
    for (std::size_t i{1}; i < std::size(zmws_); ++i) {
        const ZmwIdentity current{rgIds_[i], zmws_[i]};
        if (current == previous) {
            continue;
        }
        ordered = ordered && IdentityLess(previous, current);
        previous = current;
        ++numZmws;
    }

    if (ordered) {
        const std::uint64_t startZmw{PartitionBoundary(numZmws, chunkNum - 1, totalChunks)};
        const std::uint64_t endZmw{PartitionBoundary(numZmws, chunkNum, totalChunks)};
        if (startZmw == endZmw) {
            return {};
        }

        std::size_t startRecord{0};
        std::size_t endRecord{std::size(zmws_)};
        std::uint64_t currentZmw{0};
        previous = ZmwIdentity{rgIds_.front(), zmws_.front()};
        for (std::size_t i{1}; i < std::size(zmws_); ++i) {
            const ZmwIdentity current{rgIds_[i], zmws_[i]};
            if (current == previous) {
                continue;
            }
            previous = current;
            ++currentZmw;
            if (currentZmw == startZmw) {
                startRecord = i;
            }
            if (currentZmw == endZmw) {
                endRecord = i;
                break;
            }
        }
        return ZmwChunk{offsets_[startRecord], endRecord - startRecord, endZmw - startZmw};
    }

    // Preserve UniqueZmws()/Find() semantics for unusual indexes whose identities
    // are not ordered. One summary per unique identity replaces the two full point-
    // query hash indexes and per-record index vectors used by the old chunking path.
    struct Summary
    {
        std::int64_t FirstOffset;
        std::uint64_t NumRecords;
    };

    std::vector<Summary> summaries;
    summaries.reserve(std::size(zmws_));
    std::unordered_map<std::uint64_t, std::size_t> summaryByIdentity;
    summaryByIdentity.reserve(std::size(zmws_));
    for (std::size_t i{0}; i < std::size(zmws_); ++i) {
        const std::uint64_t key{IdentityKey(rgIds_[i], zmws_[i])};
        const auto [it, inserted]{summaryByIdentity.try_emplace(key, std::size(summaries))};
        if (inserted) {
            summaries.push_back(Summary{offsets_[i], 1});
        } else {
            ++summaries[it->second].NumRecords;
        }
    }

    numZmws = std::size(summaries);
    const std::uint64_t startZmw{PartitionBoundary(numZmws, chunkNum - 1, totalChunks)};
    const std::uint64_t endZmw{PartitionBoundary(numZmws, chunkNum, totalChunks)};
    if (startZmw == endZmw) {
        return {};
    }
    std::uint64_t numRecords{0};
    for (std::uint64_t i{startZmw}; i < endZmw; ++i) {
        numRecords += summaries[static_cast<std::size_t>(i)].NumRecords;
    }
    return ZmwChunk{summaries[static_cast<std::size_t>(startZmw)].FirstOffset, numRecords,
                    endZmw - startZmw};
}

std::int64_t ZmwIndex::FirstOffset(std::int32_t zmw) const
{
    BuildIndex();

    return FirstOffsetOrThrow(zmwIndex_, zmw, offsets_, ZmwNotFoundError{zmw});
}

std::int64_t ZmwIndex::FirstOffset(ZmwIdentity id) const
{
    BuildIndex();

    return FirstOffsetOrThrow(identityIndex_, IdentityKey(id), offsets_,
                              ZmwIdentityNotFoundError{id});
}

std::uint64_t ZmwIndex::NumRecords() const { return std::size(offsets_); }

ZmwIdentity ZmwIndex::Identity(const std::size_t row) const
{
    return ZmwIdentity{rgIds_.at(row), zmws_.at(row)};
}

std::span<const std::int32_t> ZmwIndex::QueryLengths() const { return queryLengths_; }

std::uint64_t ZmwIndex::NumZmws() const
{
    BuildIndex();
    return cachedNumZmws_;
}

}  // namespace Samoa
}  // namespace PacBio
