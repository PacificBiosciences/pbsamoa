#include <pbsamoa/index/ZmwIndex.hpp>

#include "BinaryUtils.hpp"
#include "ZmiInternal.hpp"

#include <pbsamoa/core/Bgzf.hpp>

#include <array>
#include <format>
#include <limits>
#include <optional>
#include <ranges>
#include <stdexcept>
#include <string_view>
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
        zmwIndex_ = std::move(other.zmwIndex_);
        identityIndex_ = std::move(other.identityIndex_);
        cachedNumZmws_ = other.cachedNumZmws_;
        // indexOnce_ cannot be moved; leave as default (unset).
        // call_once will rebuild if needed, but moved maps are already populated.
    }
    return *this;
}

namespace {

constexpr std::size_t ZMI_ENTRY_MIN_SIZE{16};  // rgId(4) + zmw(4) + virtualOffset(8)
constexpr std::size_t BGZF_MAX_BLOCK_SIZE{65536};

// PBI format constants
constexpr std::size_t PBI_HEADER_SIZE{32};
constexpr std::array<char, 4> PBI_MAGIC{'P', 'B', 'I', '\1'};

std::filesystem::path SidecarPath(const std::filesystem::path& bamPath, std::string_view suffix)
{
    return std::filesystem::path{bamPath.string() + std::string{suffix}};
}

std::size_t CheckedMul(std::size_t lhs, std::size_t rhs, std::string_view what)
{
    if ((lhs != 0) && (rhs > (std::numeric_limits<std::size_t>::max() / lhs))) {
        throw std::runtime_error{std::format("{} size overflows", what)};
    }
    return lhs * rhs;
}

std::size_t CheckedAdd(std::size_t lhs, std::size_t rhs, std::string_view what)
{
    if (rhs > (std::numeric_limits<std::size_t>::max() - lhs)) {
        throw std::runtime_error{std::format("{} size overflows", what)};
    }
    return lhs + rhs;
}

std::vector<std::byte> ReadAllBgzfData(BgzfReader& reader)
{
    std::vector<std::byte> data;
    std::array<std::byte, BGZF_MAX_BLOCK_SIZE> blockBuf{};

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

void ReadInt32Column(std::vector<std::int32_t>& output, std::span<const std::byte> data,
                     std::size_t offset, std::uint32_t count)
{
    output.reserve(std::size(output) + count);
    for (std::uint32_t i{0}; i < count; ++i) {
        output.push_back(ReadI32LE(std::data(data) + offset + (i * sizeof(std::int32_t))));
    }
}

void ReadInt64Column(std::vector<std::int64_t>& output, std::span<const std::byte> data,
                     std::size_t offset, std::uint32_t count)
{
    output.reserve(std::size(output) + count);
    for (std::uint32_t i{0}; i < count; ++i) {
        output.push_back(ReadI64LE(std::data(data) + offset + (i * sizeof(std::int64_t))));
    }
}

void AppendOffsets(std::vector<std::int64_t>& output, std::span<const std::ptrdiff_t> indices,
                   std::span<const std::int64_t> offsets)
{
    output.reserve(std::size(output) + std::size(indices));
    for (const std::ptrdiff_t idx : indices) {
        output.push_back(offsets[idx]);
    }
}

std::vector<std::int64_t> BuildOffsetResult(std::span<const std::ptrdiff_t> indices,
                                            std::span<const std::int64_t> offsets)
{
    std::vector<std::int64_t> result;
    AppendOffsets(result, indices, offsets);
    return result;
}

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
    if (entrySize < ZMI_ENTRY_MIN_SIZE) {
        throw std::runtime_error{std::format("ZMI entrySize too small: {}", entrySize)};
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
    }

    return index;
}

ZmwIndex ZmwIndex::FromPbi(const std::filesystem::path& path)
{
    BgzfReader reader{path};
    std::vector<std::byte> data{ReadAllBgzfData(reader)};

    // Validate minimum size for header
    if (std::size(data) < PBI_HEADER_SIZE) {
        throw std::runtime_error{"PBI file too small: missing header"};
    }

    // Validate magic "PBI\1"
    if (std::memcmp(std::data(data), std::data(PBI_MAGIC), 4) != 0) {
        throw std::runtime_error{"PBI file has invalid magic bytes"};
    }
    // Read numReads from offset 10 (uint32 LE)
    const std::uint32_t numReads{ReadU32LE(std::data(data) + 10)};

    // BasicData columns start at offset 32 (end of header).
    // Columns are stored sequentially (not interleaved):
    //   rgId:       numReads × int32
    //   qStart:     numReads × int32  (skip)
    //   qEnd:       numReads × int32  (skip)
    //   holeNumber: numReads × int32
    //   readQual:   numReads × float32 (skip)
    //   ctxtFlag:   numReads × uint8   (skip)
    //   fileOffset: numReads × int64

    const std::size_t readsAsSize{numReads};
    const std::size_t int32ColBytes{
        CheckedMul(readsAsSize, sizeof(std::int32_t), "PBI int32 column")};
    const std::size_t float32ColBytes{CheckedMul(readsAsSize, sizeof(float), "PBI float32 column")};
    const std::size_t uint8ColBytes{
        CheckedMul(readsAsSize, sizeof(std::uint8_t), "PBI uint8 column")};
    const std::size_t int64ColBytes{
        CheckedMul(readsAsSize, sizeof(std::int64_t), "PBI int64 column")};

    // Validate that we have enough data for all BasicData columns
    std::size_t requiredSize{PBI_HEADER_SIZE};
    requiredSize = CheckedAdd(requiredSize, CheckedMul(4, int32ColBytes, "PBI int32 columns"),
                              "PBI total size");
    requiredSize = CheckedAdd(requiredSize, float32ColBytes, "PBI total size");
    requiredSize = CheckedAdd(requiredSize, uint8ColBytes, "PBI total size");
    requiredSize = CheckedAdd(requiredSize, int64ColBytes, "PBI total size");
    if (std::size(data) < requiredSize) {
        throw std::runtime_error{"PBI file too small: incomplete BasicData section"};
    }

    ZmwIndex index;
    index.rgIds_.reserve(numReads);
    index.zmws_.reserve(numReads);
    index.offsets_.reserve(numReads);

    std::size_t offset{PBI_HEADER_SIZE};

    // Column 1: rgId (numReads × int32)
    ReadInt32Column(index.rgIds_, data, offset, numReads);
    offset += int32ColBytes;

    // Column 2: qStart (skip)
    offset += int32ColBytes;

    // Column 3: qEnd (skip)
    offset += int32ColBytes;

    // Column 4: holeNumber (numReads × int32)
    ReadInt32Column(index.zmws_, data, offset, numReads);
    offset += int32ColBytes;

    // Column 5: readQual (skip)
    offset += float32ColBytes;

    // Column 6: ctxtFlag (skip)
    offset += uint8ColBytes;

    // Column 7: fileOffset (numReads × int64)
    ReadInt64Column(index.offsets_, data, offset, numReads);

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

void ZmwIndex::BuildIndex() const
{
    std::call_once(indexOnce_, [this]() {
        // Skip if already populated (e.g., moved from a built index)
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
    });
}

std::vector<std::int64_t> ZmwIndex::Find(std::int32_t zmw) const
{
    BuildIndex();

    const auto it{zmwIndex_.find(zmw)};
    if (it == std::ranges::end(zmwIndex_)) {
        return {};
    }

    return BuildOffsetResult(it->second, offsets_);
}

std::vector<std::int64_t> ZmwIndex::Find(ZmwIdentity id) const
{
    BuildIndex();

    const std::uint64_t key{IdentityKey(id.rgId, id.zmw)};
    const auto it{identityIndex_.find(key)};
    if (it == std::ranges::end(identityIndex_)) {
        return {};
    }

    return BuildOffsetResult(it->second, offsets_);
}

std::vector<std::int64_t> ZmwIndex::Find(std::span<const ZmwIdentity> ids) const
{
    BuildIndex();

    std::vector<std::int64_t> result;
    result.reserve(std::size(ids));
    for (const ZmwIdentity& id : ids) {
        const std::uint64_t key{IdentityKey(id.rgId, id.zmw)};
        const auto it{identityIndex_.find(key)};
        if (it != std::ranges::end(identityIndex_)) {
            AppendOffsets(result, it->second, offsets_);
        }
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
        const std::uint64_t key{IdentityKey(id.rgId, id.zmw)};
        if (seen.insert(key).second) {
            result.push_back(id);
        }
    }

    return result;
}

std::int64_t ZmwIndex::FirstOffset(std::int32_t zmw) const
{
    BuildIndex();

    const auto it{zmwIndex_.find(zmw)};
    if (it == std::ranges::end(zmwIndex_)) {
        throw std::runtime_error{std::format("ZMW not found: {}", zmw)};
    }
    // First entry is the lowest file-order index (indices are inserted in order)
    return offsets_[it->second.front()];
}

std::int64_t ZmwIndex::FirstOffset(ZmwIdentity id) const
{
    BuildIndex();

    const std::uint64_t key{IdentityKey(id.rgId, id.zmw)};
    const auto it{identityIndex_.find(key)};
    if (it == std::ranges::end(identityIndex_)) {
        throw std::runtime_error{
            std::format("ZMW identity not found: rgId={} zmw={}", id.rgId, id.zmw)};
    }
    return offsets_[it->second.front()];
}

std::uint64_t ZmwIndex::NumRecords() const { return std::size(offsets_); }

std::uint64_t ZmwIndex::NumZmws() const
{
    BuildIndex();
    return cachedNumZmws_;
}

}  // namespace Samoa
}  // namespace PacBio
