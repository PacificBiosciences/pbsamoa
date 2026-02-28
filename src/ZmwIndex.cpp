#include <pbsamoa/index/ZmwIndex.hpp>

#include "ZmiInternal.hpp"

#include <pbsamoa/core/Bgzf.hpp>

#include <array>
#include <format>
#include <limits>
#include <optional>
#include <ranges>
#include <stdexcept>
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

template <typename T>
T ReadLE(const std::byte* src)
{
    T value{};
    std::byte* dst{reinterpret_cast<std::byte*>(&value)};
    std::ranges::copy_n(src, sizeof(value), dst);
    return value;
}

std::uint16_t ReadLE16(const std::byte* src) { return ReadLE<std::uint16_t>(src); }

std::int32_t ReadLE32Signed(const std::byte* src) { return ReadLE<std::int32_t>(src); }

std::int64_t ReadLE64Signed(const std::byte* src) { return ReadLE<std::int64_t>(src); }

std::uint32_t ReadLE32Unsigned(const std::byte* src) { return ReadLE<std::uint32_t>(src); }

}  // namespace

ZmwIndex ZmwIndex::FromZmi(const std::filesystem::path& path)
{
    BgzfReader reader{path};

    // Read all decompressed data into a single buffer.
    // BGZF blocks decompress to at most 65536 bytes each.
    std::vector<std::byte> data;
    std::array<std::byte, BGZF_MAX_BLOCK_SIZE> blockBuf{};

    while (true) {
        const std::optional<std::size_t> bytesRead{reader.ReadBlock(blockBuf)};
        if (!bytesRead.has_value() || (*bytesRead == 0)) {
            break;
        }
        data.insert(std::ranges::end(data), std::ranges::begin(blockBuf),
                    std::ranges::begin(blockBuf) + static_cast<std::ptrdiff_t>(*bytesRead));
    }

    // Validate minimum size for header
    if (std::size(data) < detail::ZMI_HEADER_SIZE) {
        throw std::runtime_error{"ZMI file too small: missing header"};
    }

    // Validate magic "ZMI\1"
    if (std::memcmp(std::data(data), std::data(detail::ZMI_MAGIC), 4) != 0) {
        throw std::runtime_error{"ZMI file has invalid magic bytes"};
    }

    // Read entrySize from offset 8 (2 bytes LE)
    const std::uint16_t entrySize{ReadLE16(std::data(data) + 8)};
    if (entrySize < ZMI_ENTRY_MIN_SIZE) {
        throw std::runtime_error{std::format("ZMI entrySize too small: {}", entrySize)};
    }

    // Iterate over entries starting at offset 64 with stride = entrySize
    const std::size_t dataSize{std::size(data)};
    const std::size_t numEntries{(dataSize - detail::ZMI_HEADER_SIZE) / entrySize};

    ZmwIndex index;
    index.rgIds_.reserve(numEntries);
    index.zmws_.reserve(numEntries);
    index.offsets_.reserve(numEntries);

    for (std::size_t i{0}; i < numEntries; ++i) {
        const std::byte* entry{std::data(data) + detail::ZMI_HEADER_SIZE + (i * entrySize)};

        const std::int32_t rgId{ReadLE32Signed(entry)};
        const std::int32_t zmw{ReadLE32Signed(entry + 4)};
        const std::int64_t virtualOffset{ReadLE64Signed(entry + 8)};

        index.rgIds_.push_back(rgId);
        index.zmws_.push_back(zmw);
        index.offsets_.push_back(virtualOffset);
    }

    return index;
}

ZmwIndex ZmwIndex::FromPbi(const std::filesystem::path& path)
{
    BgzfReader reader{path};

    // Read all decompressed data into a single buffer.
    std::vector<std::byte> data;
    std::array<std::byte, BGZF_MAX_BLOCK_SIZE> blockBuf{};

    while (true) {
        const std::optional<std::size_t> bytesRead{reader.ReadBlock(blockBuf)};
        if (!bytesRead.has_value() || (*bytesRead == 0)) {
            break;
        }
        data.insert(std::ranges::end(data), std::ranges::begin(blockBuf),
                    std::ranges::begin(blockBuf) + static_cast<std::ptrdiff_t>(*bytesRead));
    }

    // Validate minimum size for header
    if (std::size(data) < PBI_HEADER_SIZE) {
        throw std::runtime_error{"PBI file too small: missing header"};
    }

    // Validate magic "PBI\1"
    if (std::memcmp(std::data(data), std::data(PBI_MAGIC), 4) != 0) {
        throw std::runtime_error{"PBI file has invalid magic bytes"};
    }

    // Read numReads from offset 10 (uint32 LE)
    const std::uint32_t numReads{ReadLE32Unsigned(std::data(data) + 10)};

    // BasicData columns start at offset 32 (end of header).
    // Columns are stored sequentially (not interleaved):
    //   rgId:       numReads × int32
    //   qStart:     numReads × int32  (skip)
    //   qEnd:       numReads × int32  (skip)
    //   holeNumber: numReads × int32
    //   readQual:   numReads × float32 (skip)
    //   ctxtFlag:   numReads × uint8   (skip)
    //   fileOffset: numReads × int64

    const std::size_t int32ColBytes = numReads * sizeof(std::int32_t);
    const std::size_t float32ColBytes = numReads * sizeof(float);
    const std::size_t uint8ColBytes = numReads * sizeof(std::uint8_t);
    const std::size_t int64ColBytes = numReads * sizeof(std::int64_t);

    // Validate that we have enough data for all BasicData columns
    const std::size_t requiredSize{PBI_HEADER_SIZE + (4 * int32ColBytes) + float32ColBytes +
                                   uint8ColBytes + int64ColBytes};
    if (std::size(data) < requiredSize) {
        throw std::runtime_error{"PBI file too small: incomplete BasicData section"};
    }

    ZmwIndex index;
    index.rgIds_.reserve(numReads);
    index.zmws_.reserve(numReads);
    index.offsets_.reserve(numReads);

    std::size_t offset{PBI_HEADER_SIZE};

    // Column 1: rgId (numReads × int32)
    for (std::uint32_t i{0}; i < numReads; ++i) {
        index.rgIds_.push_back(
            ReadLE32Signed(std::data(data) + offset + (i * sizeof(std::int32_t))));
    }
    offset += int32ColBytes;

    // Column 2: qStart (skip)
    offset += int32ColBytes;

    // Column 3: qEnd (skip)
    offset += int32ColBytes;

    // Column 4: holeNumber (numReads × int32)
    for (std::uint32_t i{0}; i < numReads; ++i) {
        index.zmws_.push_back(
            ReadLE32Signed(std::data(data) + offset + (i * sizeof(std::int32_t))));
    }
    offset += int32ColBytes;

    // Column 5: readQual (skip)
    offset += float32ColBytes;

    // Column 6: ctxtFlag (skip)
    offset += uint8ColBytes;

    // Column 7: fileOffset (numReads × int64)
    for (std::uint32_t i{0}; i < numReads; ++i) {
        index.offsets_.push_back(
            ReadLE64Signed(std::data(data) + offset + (i * sizeof(std::int64_t))));
    }

    return index;
}

ZmwIndex ZmwIndex::Open(const std::filesystem::path& bamPath)
{
    const std::filesystem::path zmiPath{bamPath.string() + ".zmi"};
    if (std::filesystem::exists(zmiPath)) {
        return FromZmi(zmiPath);
    }
    const std::filesystem::path pbiPath{bamPath.string() + ".pbi"};
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

    std::vector<std::int64_t> result;
    result.reserve(std::size(it->second));
    for (const std::ptrdiff_t idx : it->second) {
        result.push_back(offsets_[idx]);
    }
    return result;
}

std::vector<std::int64_t> ZmwIndex::Find(ZmwIdentity id) const
{
    BuildIndex();

    const std::uint64_t key{IdentityKey(id.rgId, id.zmw)};
    const auto it{identityIndex_.find(key)};
    if (it == std::ranges::end(identityIndex_)) {
        return {};
    }

    std::vector<std::int64_t> result;
    result.reserve(std::size(it->second));
    for (const std::ptrdiff_t idx : it->second) {
        result.push_back(offsets_[idx]);
    }
    return result;
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
            for (const std::ptrdiff_t idx : it->second) {
                result.push_back(offsets_[idx]);
            }
        }
    }
    return result;
}

std::vector<ZmwIdentity> ZmwIndex::UniqueZmws() const
{
    std::vector<ZmwIdentity> result;
    result.reserve(std::size(zmws_));

    std::int32_t lastRgId{std::numeric_limits<std::int32_t>::min()};
    std::int32_t lastZmw{std::numeric_limits<std::int32_t>::min()};

    for (std::ptrdiff_t i{0}; i < std::ssize(zmws_); ++i) {
        if ((rgIds_[i] != lastRgId) || (zmws_[i] != lastZmw)) {
            lastRgId = rgIds_[i];
            lastZmw = zmws_[i];
            result.push_back(ZmwIdentity{lastRgId, lastZmw});
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
