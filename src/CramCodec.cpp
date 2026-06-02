#include <pbsamoa/cram/CramCodec.hpp>
#include <pbsamoa/cram/CramStructs.hpp>

#include <algorithm>
#include <format>
#include <stdexcept>
#include <string_view>
#include <utility>

#include <cstring>

namespace PacBio {
namespace Samoa {

namespace {

std::size_t CheckedByteArrayLength(const std::int32_t length, const std::string_view codecName)
{
    if (length < 0) {
        throw std::runtime_error{std::format("{}: negative array length {}", codecName, length)};
    }
    return static_cast<std::size_t>(length);
}

[[noreturn]] void ThrowUnsupported(const std::string_view codecName,
                                   const std::string_view operation)
{
    throw std::runtime_error{std::format("{}: {} not supported", codecName, operation)};
}

template <typename Codec>
std::byte DecodeByteViaInt(Codec& codec, CramBitReader& coreReader,
                           CramExternalBlockStore& extBlocks)
{
    return static_cast<std::byte>(codec.DecodeInt(coreReader, extBlocks));
}

template <typename Codec>
void EncodeByteViaInt(Codec& codec, const std::byte value, CramBitWriter& coreWriter,
                      CramExternalBlockStore& extBlocks)
{
    codec.EncodeInt(static_cast<std::int32_t>(static_cast<std::uint8_t>(value)), coreWriter,
                    extBlocks);
}

CramEncodingDescriptor ParseNestedEncodingDescriptor(std::span<const std::byte> parameters,
                                                     std::size_t& pos)
{
    CramEncodingDescriptor desc;
    std::size_t bytesRead = 0;
    desc.CodecId = static_cast<CramCodecId>(
        ReadItf8(std::span<const std::byte>{parameters}.subspan(pos), bytesRead));
    pos += bytesRead;

    const auto parameterLength =
        ReadItf8(std::span<const std::byte>{parameters}.subspan(pos), bytesRead);
    pos += bytesRead;

    desc.Parameters.assign(std::begin(parameters) + pos,
                           std::begin(parameters) + pos + parameterLength);
    pos += parameterLength;
    return desc;
}

}  // namespace

// ---------------------------------------------------------------------------
// CramBitReader
// ---------------------------------------------------------------------------

CramBitReader::CramBitReader(std::span<const std::byte> data) : data_{data} {}

std::int32_t CramBitReader::ReadBit()
{
    if (bytePos_ >= std::size(data_)) {
        throw std::runtime_error("CramBitReader: past end of data");
    }
    const int bit = (static_cast<std::uint8_t>(data_[bytePos_]) >> bitPos_) & 1;
    if (--bitPos_ < 0) {
        bitPos_ = 7;
        ++bytePos_;
    }
    return bit;
}

std::int32_t CramBitReader::ReadBits(int numBits)
{
    std::int32_t value = 0;
    for (int i = 0; i < numBits; ++i) {
        value = (value << 1) | ReadBit();
    }
    return value;
}

std::size_t CramBitReader::BytePosition() const { return bytePos_; }

// ---------------------------------------------------------------------------
// CramBitWriter
// ---------------------------------------------------------------------------

CramBitWriter::CramBitWriter() = default;

void CramBitWriter::WriteBit(int bit)
{
    if (bit) {
        currentByte_ |= static_cast<std::byte>(1 << bitPos_);
    }
    if (--bitPos_ < 0) {
        data_.push_back(currentByte_);
        currentByte_ = std::byte{0};
        bitPos_ = 7;
    }
}

void CramBitWriter::WriteBits(std::int32_t value, int numBits)
{
    for (int i = numBits - 1; i >= 0; --i) {
        WriteBit((value >> i) & 1);
    }
}

void CramBitWriter::Flush()
{
    if (bitPos_ < 7) {
        // Left-shift remaining bits
        data_.push_back(currentByte_);
        currentByte_ = std::byte{0};
        bitPos_ = 7;
    }
}

std::vector<std::byte> CramBitWriter::Data() const& { return data_; }

std::vector<std::byte> CramBitWriter::Data() && { return std::move(data_); }

// ---------------------------------------------------------------------------
// CramExternalBlockStore
// ---------------------------------------------------------------------------

bool CramExternalBlockStore::IsFastContentId(const std::int32_t contentId)
{
    return contentId >= 0 && contentId < static_cast<std::int32_t>(FAST_BLOCK_LIMIT);
}

std::size_t CramExternalBlockStore::FastContentIndex(const std::int32_t contentId)
{
    return static_cast<std::size_t>(contentId);
}

CramExternalBlockStore::BlockState* CramExternalBlockStore::FindBlock(std::int32_t contentId)
{
    if (BlockState* const state = FindExistingBlock(contentId)) {
        return state;
    }
    throw std::runtime_error{std::format("CramExternalBlockStore: unknown block {}", contentId)};
}

const CramExternalBlockStore::BlockState* CramExternalBlockStore::FindExistingBlock(
    std::int32_t contentId) const
{
    if (IsFastContentId(contentId)) {
        const std::size_t idx = FastContentIndex(contentId);
        if (!fastBlockUsed_[idx]) {
            return nullptr;
        }
        return &fastBlocks_[idx];
    }

    const auto it = blocks_.find(contentId);
    if (it == blocks_.end()) {
        return nullptr;
    }
    return &it->second;
}

CramExternalBlockStore::BlockState* CramExternalBlockStore::FindExistingBlock(
    std::int32_t contentId)
{
    return const_cast<BlockState*>(std::as_const(*this).FindExistingBlock(contentId));
}

CramExternalBlockStore::BlockState* CramExternalBlockStore::FindOrCreateBlock(
    std::int32_t contentId)
{
    if (IsFastContentId(contentId)) {
        const std::size_t idx = FastContentIndex(contentId);
        fastBlockUsed_[idx] = true;
        return &fastBlocks_[idx];
    }
    return &blocks_[contentId];
}

void CramExternalBlockStore::AddBlock(std::int32_t contentId, std::span<const std::byte> data)
{
    auto* state = FindOrCreateBlock(contentId);
    state->data.assign(data.begin(), data.end());
    state->readPos = 0;
}

std::span<const std::byte> CramExternalBlockStore::ReadBytesView(std::int32_t contentId,
                                                                 std::size_t n)
{
    auto* state = FindBlock(contentId);
    if (state->readPos + n > std::size(state->data)) {
        throw std::runtime_error{
            std::format("CramExternalBlockStore: read past end of block {}", contentId)};
    }
    const auto view = std::span<const std::byte>{state->data}.subspan(state->readPos, n);
    state->readPos += n;
    return view;
}

std::vector<std::byte> CramExternalBlockStore::ReadBytes(std::int32_t contentId, std::size_t n)
{
    const auto view = ReadBytesView(contentId, n);
    return {view.begin(), view.end()};
}

std::span<const std::byte> CramExternalBlockStore::ReadBytesUntilStop(std::int32_t contentId,
                                                                      std::byte stopByte)
{
    auto* state = FindBlock(contentId);
    const auto* start = state->data.data() + state->readPos;
    const auto remaining = std::size(state->data) - state->readPos;
    const auto* found = static_cast<const std::byte*>(
        std::memchr(start, std::to_integer<int>(stopByte), remaining));
    if (!found) {
        throw std::runtime_error{
            std::format("CramExternalBlockStore: stop byte not found in block {}", contentId)};
    }
    const auto len = static_cast<std::size_t>(found - start);
    const auto view = std::span<const std::byte>{start, len};
    state->readPos += len + 1;  // advance past stop byte
    return view;
}

std::byte CramExternalBlockStore::ReadByte(std::int32_t contentId)
{
    auto* state = FindBlock(contentId);
    if (state->readPos >= std::size(state->data)) {
        throw std::runtime_error{
            std::format("CramExternalBlockStore: read past end of block {}", contentId)};
    }
    const std::size_t readPos{state->readPos};
    ++state->readPos;
    return state->data[readPos];
}

std::int32_t CramExternalBlockStore::ReadItf8(std::int32_t contentId)
{
    auto* state = FindBlock(contentId);
    std::size_t bytesRead = 0;
    const auto value = PacBio::Samoa::ReadItf8(
        std::span<const std::byte>{state->data}.subspan(state->readPos), bytesRead);
    state->readPos += bytesRead;
    return value;
}

void CramExternalBlockStore::WriteBytes(std::int32_t contentId, std::span<const std::byte> data)
{
    auto* state = FindOrCreateBlock(contentId);
    state->data.insert(state->data.end(), data.begin(), data.end());
}

void CramExternalBlockStore::WriteByte(std::int32_t contentId, std::byte value)
{
    auto* state = FindOrCreateBlock(contentId);
    state->data.push_back(value);
}

void CramExternalBlockStore::WriteItf8(std::int32_t contentId, std::int32_t value)
{
    auto* state = FindOrCreateBlock(contentId);
    PacBio::Samoa::WriteItf8(state->data, value);
}

void CramExternalBlockStore::ReserveBlock(std::int32_t contentId, std::size_t capacity)
{
    auto* state = FindOrCreateBlock(contentId);
    state->data.reserve(capacity);
}

std::span<const std::byte> CramExternalBlockStore::GetBlockData(std::int32_t contentId) const
{
    const BlockState* const state = FindExistingBlock(contentId);
    if (!state) {
        return {};
    }
    return state->data;
}

std::vector<std::byte> CramExternalBlockStore::TakeBlockData(std::int32_t contentId)
{
    BlockState* const state = FindExistingBlock(contentId);
    if (!state) {
        return {};
    }
    return std::move(state->data);
}

std::vector<std::int32_t> CramExternalBlockStore::ContentIds() const
{
    std::vector<std::int32_t> ids;
    ids.reserve(std::size(blocks_) + FAST_BLOCK_LIMIT);
    for (std::size_t i = 0; i < FAST_BLOCK_LIMIT; ++i) {
        if (fastBlockUsed_[i]) {
            ids.push_back(static_cast<std::int32_t>(i));
        }
    }
    for (const auto& [id, _] : blocks_) {
        ids.push_back(id);
    }
    std::ranges::sort(ids);
    return ids;
}

void CramExternalBlockStore::ResetPositions()
{
    for (std::size_t i = 0; i < FAST_BLOCK_LIMIT; ++i) {
        if (fastBlockUsed_[i]) {
            fastBlocks_[i].readPos = 0;
        }
    }
    for (auto& [_, state] : blocks_) {
        state.readPos = 0;
    }
}

// ---------------------------------------------------------------------------
// CramCodec — default DecodeByteArrayView (copies)
// ---------------------------------------------------------------------------

std::span<const std::byte> CramCodec::DecodeByteArrayView(CramBitReader& coreReader,
                                                          CramExternalBlockStore& extBlocks)
{
    viewScratch_ = DecodeByteArray(coreReader, extBlocks);
    return viewScratch_;
}

// ---------------------------------------------------------------------------
// NullCodec
// ---------------------------------------------------------------------------

CramCodecDecodeKind NullCodec::DecodeKind() const { return CramCodecDecodeKind::SCALAR; }

std::int32_t NullCodec::DecodeInt(CramBitReader& /*coreReader*/,
                                  CramExternalBlockStore& /*extBlocks*/)
{
    return 0;
}

std::byte NullCodec::DecodeByte(CramBitReader& /*coreReader*/,
                                CramExternalBlockStore& /*extBlocks*/)
{
    return std::byte{0};
}

std::vector<std::byte> NullCodec::DecodeByteArray(CramBitReader& /*coreReader*/,
                                                  CramExternalBlockStore& /*extBlocks*/)
{
    return {};
}

void NullCodec::EncodeInt(std::int32_t /*value*/, CramBitWriter& /*coreWriter*/,
                          CramExternalBlockStore& /*extBlocks*/)
{
}

void NullCodec::EncodeByte(std::byte /*value*/, CramBitWriter& /*coreWriter*/,
                           CramExternalBlockStore& /*extBlocks*/)
{
}

void NullCodec::EncodeByteArray(std::span<const std::byte> /*data*/, CramBitWriter& /*coreWriter*/,
                                CramExternalBlockStore& /*extBlocks*/)
{
}

// ---------------------------------------------------------------------------
// ExternalCodec
// ---------------------------------------------------------------------------

ExternalCodec::ExternalCodec(std::int32_t blockContentId) : blockContentId_{blockContentId} {}

CramCodecDecodeKind ExternalCodec::DecodeKind() const { return CramCodecDecodeKind::SCALAR; }

std::int32_t ExternalCodec::DecodeInt(CramBitReader& /*coreReader*/,
                                      CramExternalBlockStore& extBlocks)
{
    return extBlocks.ReadItf8(blockContentId_);
}

std::byte ExternalCodec::DecodeByte(CramBitReader& /*coreReader*/,
                                    CramExternalBlockStore& extBlocks)
{
    return extBlocks.ReadByte(blockContentId_);
}

std::vector<std::byte> ExternalCodec::DecodeByteArray(CramBitReader& /*coreReader*/,
                                                      CramExternalBlockStore& extBlocks)
{
    const std::int32_t len{extBlocks.ReadItf8(blockContentId_)};
    const std::size_t byteCount{CheckedByteArrayLength(len, "ExternalCodec")};
    const auto view = extBlocks.ReadBytesView(blockContentId_, byteCount);
    return {view.begin(), view.end()};
}

void ExternalCodec::EncodeInt(std::int32_t value, CramBitWriter& /*coreWriter*/,
                              CramExternalBlockStore& extBlocks)
{
    extBlocks.WriteItf8(blockContentId_, value);
}

void ExternalCodec::EncodeByte(std::byte value, CramBitWriter& /*coreWriter*/,
                               CramExternalBlockStore& extBlocks)
{
    extBlocks.WriteByte(blockContentId_, value);
}

void ExternalCodec::EncodeByteArray(std::span<const std::byte> data, CramBitWriter& /*coreWriter*/,
                                    CramExternalBlockStore& extBlocks)
{
    extBlocks.WriteItf8(blockContentId_, static_cast<std::int32_t>(std::size(data)));
    extBlocks.WriteBytes(blockContentId_, data);
}

std::int32_t ExternalCodec::BlockContentId() const { return blockContentId_; }

// ---------------------------------------------------------------------------
// HuffmanCodec
// ---------------------------------------------------------------------------

HuffmanCodec::HuffmanCodec(std::vector<std::int32_t> symbols, std::vector<std::int32_t> bitLengths)
{
    if (std::size(symbols) != std::size(bitLengths)) {
        throw std::runtime_error("HuffmanCodec: symbols/bitLengths size mismatch");
    }
    entries_.reserve(std::size(symbols));
    for (std::size_t i = 0; i < std::size(symbols); ++i) {
        entries_.push_back({.symbol = symbols[i], .code = 0, .bitLength = bitLengths[i]});
    }
    BuildCodes();
}

CramCodecDecodeKind HuffmanCodec::DecodeKind() const { return CramCodecDecodeKind::SCALAR; }

void HuffmanCodec::BuildCodes()
{
    // Sort by bit length then by symbol
    std::ranges::sort(entries_, {}, &HuffmanEntry::symbol);
    std::ranges::stable_sort(entries_, {}, &HuffmanEntry::bitLength);

    for (const auto& e : entries_) {
        if (e.bitLength < 0 || e.bitLength > 31) {
            throw std::runtime_error("HuffmanCodec: invalid bit length");
        }
    }

    decodeBuckets_.clear();
    maxBitLength_ = 0;

    // Canonical Huffman code assignment.
    std::uint32_t code = 0;
    std::int32_t prevLen = 0;
    for (auto& e : entries_) {
        if (e.bitLength > 0) {
            code <<= (e.bitLength - prevLen);
            e.code = code;
            ++code;
            prevLen = e.bitLength;
            maxBitLength_ = std::max(maxBitLength_, e.bitLength);
        }
    }

    if (maxBitLength_ == 0) {
        return;
    }

    decodeBuckets_.assign(static_cast<std::size_t>(maxBitLength_ + 1), {});
    for (std::size_t i = 0; i < std::size(entries_); ++i) {
        const auto& e = entries_[i];
        if (e.bitLength <= 0) {
            continue;
        }

        auto& bucket = decodeBuckets_[static_cast<std::size_t>(e.bitLength)];
        if (!bucket.hasEntries) {
            bucket.hasEntries = true;
            bucket.firstIndex = i;
            bucket.minCode = e.code;
            bucket.maxCode = e.code;
        } else {
            bucket.maxCode = e.code;
        }
    }
}

std::int32_t HuffmanCodec::DecodeInt(CramBitReader& coreReader,
                                     CramExternalBlockStore& /*extBlocks*/)
{
    // Single-symbol shortcut
    if (std::size(entries_) == 1) {
        return entries_[0].symbol;
    }

    std::uint32_t code = 0;
    for (std::int32_t bits = 1; bits <= maxBitLength_; ++bits) {
        code = (code << 1) | static_cast<std::uint32_t>(coreReader.ReadBit());
        const auto& bucket = decodeBuckets_[static_cast<std::size_t>(bits)];
        if (!bucket.hasEntries || code < bucket.minCode || code > bucket.maxCode) {
            continue;
        }

        const std::size_t index =
            bucket.firstIndex + static_cast<std::size_t>(code - bucket.minCode);
        if (index < std::size(entries_)) {
            const auto& entry = entries_[index];
            if (entry.bitLength == bits && entry.code == code) {
                return entry.symbol;
            }
        }
    }

    throw std::runtime_error("HuffmanCodec: decode failed — no matching code");
}

std::byte HuffmanCodec::DecodeByte(CramBitReader& coreReader, CramExternalBlockStore& extBlocks)
{
    return DecodeByteViaInt(*this, coreReader, extBlocks);
}

std::vector<std::byte> HuffmanCodec::DecodeByteArray(CramBitReader& /*coreReader*/,
                                                     CramExternalBlockStore& /*extBlocks*/)
{
    ThrowUnsupported("HuffmanCodec", "DecodeByteArray");
}

void HuffmanCodec::EncodeInt(std::int32_t value, CramBitWriter& coreWriter,
                             CramExternalBlockStore& /*extBlocks*/)
{
    for (const auto& e : entries_) {
        if (e.symbol == value) {
            if (e.bitLength > 0) {
                coreWriter.WriteBits(e.code, e.bitLength);
            }
            return;
        }
    }
    throw std::runtime_error("HuffmanCodec: symbol not in table");
}

void HuffmanCodec::EncodeByte(std::byte value, CramBitWriter& coreWriter,
                              CramExternalBlockStore& extBlocks)
{
    EncodeByteViaInt(*this, value, coreWriter, extBlocks);
}

void HuffmanCodec::EncodeByteArray(std::span<const std::byte> /*data*/,
                                   CramBitWriter& /*coreWriter*/,
                                   CramExternalBlockStore& /*extBlocks*/)
{
    ThrowUnsupported("HuffmanCodec", "EncodeByteArray");
}

// ---------------------------------------------------------------------------
// BetaCodec
// ---------------------------------------------------------------------------

BetaCodec::BetaCodec(std::int32_t offset, std::int32_t numBits) : offset_{offset}, numBits_{numBits}
{
}

CramCodecDecodeKind BetaCodec::DecodeKind() const { return CramCodecDecodeKind::SCALAR; }

std::int32_t BetaCodec::DecodeInt(CramBitReader& coreReader, CramExternalBlockStore& /*extBlocks*/)
{
    return coreReader.ReadBits(numBits_) - offset_;
}

std::byte BetaCodec::DecodeByte(CramBitReader& coreReader, CramExternalBlockStore& extBlocks)
{
    return DecodeByteViaInt(*this, coreReader, extBlocks);
}

std::vector<std::byte> BetaCodec::DecodeByteArray(CramBitReader& /*coreReader*/,
                                                  CramExternalBlockStore& /*extBlocks*/)
{
    ThrowUnsupported("BetaCodec", "DecodeByteArray");
}

void BetaCodec::EncodeInt(std::int32_t value, CramBitWriter& coreWriter,
                          CramExternalBlockStore& /*extBlocks*/)
{
    coreWriter.WriteBits(value + offset_, numBits_);
}

void BetaCodec::EncodeByte(std::byte value, CramBitWriter& coreWriter,
                           CramExternalBlockStore& extBlocks)
{
    EncodeByteViaInt(*this, value, coreWriter, extBlocks);
}

void BetaCodec::EncodeByteArray(std::span<const std::byte> /*data*/, CramBitWriter& /*coreWriter*/,
                                CramExternalBlockStore& /*extBlocks*/)
{
    ThrowUnsupported("BetaCodec", "EncodeByteArray");
}

// ---------------------------------------------------------------------------
// GammaCodec
// ---------------------------------------------------------------------------

GammaCodec::GammaCodec(std::int32_t offset) : offset_{offset} {}

CramCodecDecodeKind GammaCodec::DecodeKind() const { return CramCodecDecodeKind::SCALAR; }

std::int32_t GammaCodec::DecodeInt(CramBitReader& coreReader, CramExternalBlockStore& /*extBlocks*/)
{
    // Count leading zeros
    int numZeros = 0;
    while (coreReader.ReadBit() == 0) {
        ++numZeros;
    }
    // Read numZeros more bits
    std::int32_t value = 1;
    for (int i = 0; i < numZeros; ++i) {
        value = (value << 1) | coreReader.ReadBit();
    }
    return value - offset_;
}

std::byte GammaCodec::DecodeByte(CramBitReader& coreReader, CramExternalBlockStore& extBlocks)
{
    return DecodeByteViaInt(*this, coreReader, extBlocks);
}

std::vector<std::byte> GammaCodec::DecodeByteArray(CramBitReader& /*coreReader*/,
                                                   CramExternalBlockStore& /*extBlocks*/)
{
    ThrowUnsupported("GammaCodec", "DecodeByteArray");
}

void GammaCodec::EncodeInt(std::int32_t value, CramBitWriter& coreWriter,
                           CramExternalBlockStore& /*extBlocks*/)
{
    const auto adjusted = value + offset_;
    if (adjusted <= 0) {
        throw std::runtime_error("GammaCodec: value + offset must be > 0");
    }

    // Count bits needed
    int bits = 0;
    auto tmp = adjusted;
    while (tmp > 0) {
        ++bits;
        tmp >>= 1;
    }

    // Write (bits-1) zeros
    for (int i = 0; i < bits - 1; ++i) {
        coreWriter.WriteBit(0);
    }
    // Write the value itself in 'bits' bits
    coreWriter.WriteBits(adjusted, bits);
}

void GammaCodec::EncodeByte(std::byte value, CramBitWriter& coreWriter,
                            CramExternalBlockStore& extBlocks)
{
    EncodeByteViaInt(*this, value, coreWriter, extBlocks);
}

void GammaCodec::EncodeByteArray(std::span<const std::byte> /*data*/, CramBitWriter& /*coreWriter*/,
                                 CramExternalBlockStore& /*extBlocks*/)
{
    ThrowUnsupported("GammaCodec", "EncodeByteArray");
}

// ---------------------------------------------------------------------------
// SubexpCodec
// ---------------------------------------------------------------------------

SubexpCodec::SubexpCodec(std::int32_t offset, std::int32_t k) : offset_{offset}, k_{k} {}

CramCodecDecodeKind SubexpCodec::DecodeKind() const { return CramCodecDecodeKind::SCALAR; }

std::int32_t SubexpCodec::DecodeInt(CramBitReader& coreReader,
                                    CramExternalBlockStore& /*extBlocks*/)
{
    // Read the unary prefix: a run of 1s terminated by a 0 (the 0 is consumed here).
    std::int32_t unary = 0;
    while (coreReader.ReadBit() == 1) {
        ++unary;
    }

    // Canonical Golomb-subexponential code (htslib cram_subexp_decode, cram_codecs.c:2434):
    // an empty prefix has a k-bit suffix; a prefix of length u has a (k+u-1)-bit suffix plus
    // an added base of 2^(k+u-1).
    std::int32_t value = 0;
    if (unary == 0) {
        if (k_ > 0) {
            value = coreReader.ReadBits(k_);
        }
    } else {
        const int numBits = k_ + unary - 1;
        if (numBits > 0) {
            value = coreReader.ReadBits(numBits);
        }
        value += (std::int32_t{1} << (k_ + unary - 1));
    }

    return value - offset_;
}

std::byte SubexpCodec::DecodeByte(CramBitReader& coreReader, CramExternalBlockStore& extBlocks)
{
    return DecodeByteViaInt(*this, coreReader, extBlocks);
}

std::vector<std::byte> SubexpCodec::DecodeByteArray(CramBitReader& /*coreReader*/,
                                                    CramExternalBlockStore& /*extBlocks*/)
{
    ThrowUnsupported("SubexpCodec", "DecodeByteArray");
}

void SubexpCodec::EncodeInt(std::int32_t value, CramBitWriter& coreWriter,
                            CramExternalBlockStore& /*extBlocks*/)
{
    const std::int32_t adjusted = value + offset_;

    // Inverse of DecodeInt. The prefix length u is the smallest count with adjusted < 2^(k+u);
    // u == 0 maps to a k-bit suffix, otherwise a (k+u-1)-bit suffix over base 2^(k+u-1).
    std::int32_t unary = 0;
    while (adjusted >= (std::int64_t{1} << (k_ + unary))) {
        ++unary;
    }

    // Write the unary prefix (u ones followed by a terminating 0).
    for (std::int32_t i = 0; i < unary; ++i) {
        coreWriter.WriteBit(1);
    }
    coreWriter.WriteBit(0);

    // Write the remainder.
    if (unary == 0) {
        if (k_ > 0) {
            coreWriter.WriteBits(adjusted, k_);
        }
    } else {
        const int numBits = k_ + unary - 1;
        if (numBits > 0) {
            coreWriter.WriteBits(adjusted - (std::int32_t{1} << (k_ + unary - 1)), numBits);
        }
    }
}

void SubexpCodec::EncodeByte(std::byte value, CramBitWriter& coreWriter,
                             CramExternalBlockStore& extBlocks)
{
    EncodeByteViaInt(*this, value, coreWriter, extBlocks);
}

void SubexpCodec::EncodeByteArray(std::span<const std::byte> /*data*/,
                                  CramBitWriter& /*coreWriter*/,
                                  CramExternalBlockStore& /*extBlocks*/)
{
    ThrowUnsupported("SubexpCodec", "EncodeByteArray");
}

// ---------------------------------------------------------------------------
// ByteArrayLenCodec
// ---------------------------------------------------------------------------

ByteArrayLenCodec::ByteArrayLenCodec(std::unique_ptr<CramCodec> lenCodec,
                                     std::unique_ptr<CramCodec> dataCodec)
    : lenCodec_{std::move(lenCodec)}, dataCodec_{std::move(dataCodec)}
{
}

CramCodecDecodeKind ByteArrayLenCodec::DecodeKind() const
{
    return CramCodecDecodeKind::BYTE_ARRAY;
}

std::int32_t ByteArrayLenCodec::DecodeInt(CramBitReader& /*coreReader*/,
                                          CramExternalBlockStore& /*extBlocks*/)
{
    ThrowUnsupported("ByteArrayLenCodec", "DecodeInt");
}

std::byte ByteArrayLenCodec::DecodeByte(CramBitReader& /*coreReader*/,
                                        CramExternalBlockStore& /*extBlocks*/)
{
    ThrowUnsupported("ByteArrayLenCodec", "DecodeByte");
}

std::vector<std::byte> ByteArrayLenCodec::DecodeByteArray(CramBitReader& coreReader,
                                                          CramExternalBlockStore& extBlocks)
{
    const std::int32_t len{lenCodec_->DecodeInt(coreReader, extBlocks)};
    const std::size_t byteCount{CheckedByteArrayLength(len, "ByteArrayLenCodec")};

    // Fast path: bulk read when data codec is External (avoids per-byte virtual
    // dispatch)
    if (const auto* ext = dynamic_cast<const ExternalCodec*>(dataCodec_.get())) {
        const auto view = extBlocks.ReadBytesView(ext->BlockContentId(), byteCount);
        return {view.begin(), view.end()};
    }

    std::vector<std::byte> result(byteCount);
    for (std::int32_t i{0}; i < len; ++i) {
        result[i] = dataCodec_->DecodeByte(coreReader, extBlocks);
    }
    return result;
}

std::span<const std::byte> ByteArrayLenCodec::DecodeByteArrayView(CramBitReader& coreReader,
                                                                  CramExternalBlockStore& extBlocks)
{
    const std::int32_t len{lenCodec_->DecodeInt(coreReader, extBlocks)};
    const std::size_t byteCount{CheckedByteArrayLength(len, "ByteArrayLenCodec")};

    // Fast path: zero-copy view when data codec is External
    if (const auto* ext = dynamic_cast<const ExternalCodec*>(dataCodec_.get())) {
        return extBlocks.ReadBytesView(ext->BlockContentId(), byteCount);
    }

    // Slow path: per-byte decode into scratch buffer
    viewScratch2_.resize(byteCount);
    for (std::int32_t i{0}; i < len; ++i) {
        viewScratch2_[i] = dataCodec_->DecodeByte(coreReader, extBlocks);
    }
    return viewScratch2_;
}

void ByteArrayLenCodec::EncodeInt(std::int32_t /*value*/, CramBitWriter& /*coreWriter*/,
                                  CramExternalBlockStore& /*extBlocks*/)
{
    ThrowUnsupported("ByteArrayLenCodec", "EncodeInt");
}

void ByteArrayLenCodec::EncodeByte(std::byte /*value*/, CramBitWriter& /*coreWriter*/,
                                   CramExternalBlockStore& /*extBlocks*/)
{
    ThrowUnsupported("ByteArrayLenCodec", "EncodeByte");
}

void ByteArrayLenCodec::EncodeByteArray(std::span<const std::byte> data, CramBitWriter& coreWriter,
                                        CramExternalBlockStore& extBlocks)
{
    lenCodec_->EncodeInt(static_cast<std::int32_t>(std::size(data)), coreWriter, extBlocks);

    // Fast path: bulk write when data codec is External (avoids per-byte virtual
    // dispatch)
    if (const auto* ext = dynamic_cast<const ExternalCodec*>(dataCodec_.get())) {
        extBlocks.WriteBytes(ext->BlockContentId(), data);
        return;
    }

    for (const auto b : data) {
        dataCodec_->EncodeByte(b, coreWriter, extBlocks);
    }
}

// ---------------------------------------------------------------------------
// ByteArrayStopCodec
// ---------------------------------------------------------------------------

ByteArrayStopCodec::ByteArrayStopCodec(std::byte stopByte, std::int32_t blockContentId)
    : stopByte_{stopByte}, blockContentId_{blockContentId}
{
}

CramCodecDecodeKind ByteArrayStopCodec::DecodeKind() const
{
    return CramCodecDecodeKind::BYTE_ARRAY;
}

std::int32_t ByteArrayStopCodec::DecodeInt(CramBitReader& /*coreReader*/,
                                           CramExternalBlockStore& /*extBlocks*/)
{
    ThrowUnsupported("ByteArrayStopCodec", "DecodeInt");
}

std::byte ByteArrayStopCodec::DecodeByte(CramBitReader& /*coreReader*/,
                                         CramExternalBlockStore& /*extBlocks*/)
{
    ThrowUnsupported("ByteArrayStopCodec", "DecodeByte");
}

std::vector<std::byte> ByteArrayStopCodec::DecodeByteArray(CramBitReader& /*coreReader*/,
                                                           CramExternalBlockStore& extBlocks)
{
    const auto view = extBlocks.ReadBytesUntilStop(blockContentId_, stopByte_);
    return {view.begin(), view.end()};
}

std::span<const std::byte> ByteArrayStopCodec::DecodeByteArrayView(
    CramBitReader& /*coreReader*/, CramExternalBlockStore& extBlocks)
{
    return extBlocks.ReadBytesUntilStop(blockContentId_, stopByte_);
}

void ByteArrayStopCodec::EncodeInt(std::int32_t /*value*/, CramBitWriter& /*coreWriter*/,
                                   CramExternalBlockStore& /*extBlocks*/)
{
    ThrowUnsupported("ByteArrayStopCodec", "EncodeInt");
}

void ByteArrayStopCodec::EncodeByte(std::byte /*value*/, CramBitWriter& /*coreWriter*/,
                                    CramExternalBlockStore& /*extBlocks*/)
{
    ThrowUnsupported("ByteArrayStopCodec", "EncodeByte");
}

void ByteArrayStopCodec::EncodeByteArray(std::span<const std::byte> data,
                                         CramBitWriter& /*coreWriter*/,
                                         CramExternalBlockStore& extBlocks)
{
    extBlocks.WriteBytes(blockContentId_, data);
    const std::array<std::byte, 1> arr{stopByte_};
    extBlocks.WriteBytes(blockContentId_, arr);
}

// ---------------------------------------------------------------------------
// Codec factory
// ---------------------------------------------------------------------------

std::unique_ptr<CramCodec> CreateCodec(const CramEncodingDescriptor& desc)
{
    switch (desc.CodecId) {
        case CramCodecId::NONE:
            return std::make_unique<NullCodec>();

        case CramCodecId::EXTERNAL: {
            std::size_t n = 0;
            const auto blockId = ReadItf8(desc.Parameters, n);
            return std::make_unique<ExternalCodec>(blockId);
        }

        case CramCodecId::HUFFMAN: {
            std::size_t pos = 0;
            std::size_t n = 0;

            // Read symbols array
            const auto numSymbols =
                ReadItf8(std::span<const std::byte>{desc.Parameters}.subspan(pos), n);
            pos += n;
            std::vector<std::int32_t> symbols(numSymbols);
            for (std::int32_t i = 0; i < numSymbols; ++i) {
                symbols[i] = ReadItf8(std::span<const std::byte>{desc.Parameters}.subspan(pos), n);
                pos += n;
            }

            // Read bit lengths array
            const auto numLengths =
                ReadItf8(std::span<const std::byte>{desc.Parameters}.subspan(pos), n);
            pos += n;
            std::vector<std::int32_t> bitLengths(numLengths);
            for (std::int32_t i = 0; i < numLengths; ++i) {
                bitLengths[i] =
                    ReadItf8(std::span<const std::byte>{desc.Parameters}.subspan(pos), n);
                pos += n;
            }

            return std::make_unique<HuffmanCodec>(std::move(symbols), std::move(bitLengths));
        }

        case CramCodecId::BYTE_ARRAY_LEN: {
            std::size_t pos = 0;
            auto lenDesc = ParseNestedEncodingDescriptor(desc.Parameters, pos);
            auto dataDesc = ParseNestedEncodingDescriptor(desc.Parameters, pos);

            return std::make_unique<ByteArrayLenCodec>(CreateCodec(lenDesc), CreateCodec(dataDesc));
        }

        case CramCodecId::BYTE_ARRAY_STOP: {
            if (std::size(desc.Parameters) < 2) {
                throw std::runtime_error("ByteArrayStop: insufficient parameters");
            }
            const auto stopByte = desc.Parameters[0];
            std::size_t n = 0;
            const auto blockId =
                ReadItf8(std::span<const std::byte>{desc.Parameters}.subspan(1), n);
            return std::make_unique<ByteArrayStopCodec>(stopByte, blockId);
        }

        case CramCodecId::BETA: {
            std::size_t pos = 0;
            std::size_t n = 0;
            const auto offset =
                ReadItf8(std::span<const std::byte>{desc.Parameters}.subspan(pos), n);
            pos += n;
            const auto numBits =
                ReadItf8(std::span<const std::byte>{desc.Parameters}.subspan(pos), n);
            return std::make_unique<BetaCodec>(offset, numBits);
        }

        case CramCodecId::SUBEXP: {
            std::size_t pos = 0;
            std::size_t n = 0;
            const auto offset =
                ReadItf8(std::span<const std::byte>{desc.Parameters}.subspan(pos), n);
            pos += n;
            const auto k = ReadItf8(std::span<const std::byte>{desc.Parameters}.subspan(pos), n);
            return std::make_unique<SubexpCodec>(offset, k);
        }

        case CramCodecId::GAMMA: {
            std::size_t n = 0;
            const auto offset = ReadItf8(desc.Parameters, n);
            return std::make_unique<GammaCodec>(offset);
        }

        default:
            throw std::runtime_error{
                std::format("CreateCodec: unknown codec ID {}", std::to_underlying(desc.CodecId))};
    }
}

}  // namespace Samoa
}  // namespace PacBio
