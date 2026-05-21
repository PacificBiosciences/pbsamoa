#ifndef PBSAMOA_TOOLS_ZMI_INDEX_RECORDPARSER_HPP
#define PBSAMOA_TOOLS_ZMI_INDEX_RECORDPARSER_HPP

#include <pbsamoa/core/Bgzf.hpp>
#include <pbsamoa/core/Endian.hpp>
#include <pbsamoa/core/Tags.hpp>

#include "../../ZmwUtils.hpp"

#include <algorithm>
#include <deque>
#include <optional>
#include <span>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

#include <cstddef>
#include <cstdint>
#include <cstring>

namespace PacBio {
namespace Samoa {
namespace ZmiIndex {

struct RecordEntry
{
    std::int32_t rgId;
    std::int32_t zmw;
    VirtualOffset virtualOffset;
};

namespace detail {

/// \brief Compute length in bytes of the value portion of a BAM aux field by
/// type code. \p data must point at the first byte of the value (after the
/// 2-byte tag key and 1-byte type code). Returns the consumed value length, or
/// std::nullopt on unknown type / truncated buffer.
inline std::optional<std::size_t> AuxValueLength(std::byte typeByte,
                                                 std::span<const std::byte> data)
{
    const char type{static_cast<char>(typeByte)};
    switch (type) {
        case 'A':
        case 'c':
        case 'C':
            return 1U;
        case 's':
        case 'S':
            return 2U;
        case 'i':
        case 'I':
        case 'f':
            return 4U;
        case 'd':
            return 8U;
        case 'Z':
        case 'H': {
            for (std::size_t i{0}; i < std::size(data); ++i) {
                if (data[i] == std::byte{0}) {
                    return i + 1U;
                }
            }
            return std::nullopt;
        }
        case 'B': {
            if (std::size(data) < 5U) {
                return std::nullopt;
            }
            const char subtype{static_cast<char>(data[0])};
            const std::uint32_t count{ReadU32LE(std::data(data) + 1U)};
            std::size_t elemSize{0};
            switch (subtype) {
                case 'c':
                case 'C':
                    elemSize = 1U;
                    break;
                case 's':
                case 'S':
                    elemSize = 2U;
                    break;
                case 'i':
                case 'I':
                case 'f':
                    elemSize = 4U;
                    break;
                default:
                    return std::nullopt;
            }
            return 5U + (static_cast<std::size_t>(count) * elemSize);
        }
        default:
            return std::nullopt;
    }
}

/// \brief Walk a BAM aux block and return the value of the \p target tag as a
/// string_view (only meaningful when target is a Z-typed tag, e.g., RG). Returns
/// nullopt if not found / malformed.
inline std::optional<std::string_view> FindZTag(std::span<const std::byte> aux, TagKey target)
{
    std::size_t i{0};
    while ((i + 3U) <= std::size(aux)) {
        const TagKey key{static_cast<char>(aux[i]), static_cast<char>(aux[i + 1U])};
        const std::byte typeByte{aux[i + 2U]};
        const std::span<const std::byte> rest{aux.subspan(i + 3U)};
        const std::optional<std::size_t> valueLen{AuxValueLength(typeByte, rest)};
        if (!valueLen) {
            return std::nullopt;
        }
        if (key == target) {
            if (static_cast<char>(typeByte) != 'Z') {
                return std::nullopt;
            }
            if (*valueLen == 0U) {
                return std::string_view{};
            }
            return std::string_view{reinterpret_cast<const char*>(std::data(rest)), *valueLen - 1U};
        }
        i += 3U + *valueLen;
    }
    return std::nullopt;
}

/// \brief Parse (rgId, zmw) for a BAM record body. \p recordBody is the bytes
/// following the 4-byte block_size prefix; its size equals block_size value.
inline std::pair<std::int32_t, std::int32_t> ParseRecordIdentity(
    std::span<const std::byte> recordBody)
{
    constexpr std::size_t FIXED_FIELDS_SIZE{32U};
    if (std::size(recordBody) < FIXED_FIELDS_SIZE) {
        throw std::runtime_error{"zmi-index: BAM record smaller than fixed header"};
    }

    const std::uint8_t nameLen{static_cast<std::uint8_t>(recordBody[8])};
    const std::uint16_t cigarOpCount{ReadU16LE(std::data(recordBody) + 12U)};
    const std::uint32_t seqLength{ReadU32LE(std::data(recordBody) + 16U)};

    if (nameLen == 0U) {
        throw std::runtime_error{"zmi-index: BAM record has zero l_read_name"};
    }

    const std::size_t nameOffset{FIXED_FIELDS_SIZE};
    const std::size_t cigarOffset{nameOffset + nameLen};
    const std::size_t seqOffset{cigarOffset + (std::size_t{4U} * cigarOpCount)};
    const std::size_t qualOffset{seqOffset + ((seqLength + 1U) / 2U)};
    const std::size_t auxOffset{qualOffset + seqLength};

    if (auxOffset > std::size(recordBody)) {
        throw std::runtime_error{"zmi-index: BAM record variable-length fields exceed block_size"};
    }
    if (recordBody[nameOffset + nameLen - 1U] != std::byte{0}) {
        throw std::runtime_error{"zmi-index: BAM record name not NUL-terminated"};
    }

    const std::string_view name{reinterpret_cast<const char*>(std::data(recordBody) + nameOffset),
                                static_cast<std::size_t>(nameLen) - 1U};
    const std::span<const std::byte> aux{recordBody.subspan(auxOffset)};

    const std::optional<std::string_view> rgText{FindZTag(aux, RG_TAG)};
    const std::int32_t rgId{rgText ? ParseReadGroupId(*rgText) : 0};
    const std::int32_t zmw{ParseZmwFromName(name)};
    return {rgId, zmw};
}

}  // namespace detail

/// \brief Streams decompressed BGZF bytes (BAM header + records) and emits
/// per-record `(rgId, zmw, virtualOffset)` index entries.
class RecordParser
{
public:
    /// \brief Append a decompressed BGZF block to the parser's buffer.
    /// \param[in] fileOffset original file position of this block (used for VO)
    /// \param[in] data decompressed bytes
    void Feed(std::uint64_t fileOffset, std::span<const std::byte> data);

    /// \brief Attempt to emit the next record. Returns nullopt if more data is
    /// needed. Throws std::runtime_error on malformed input.
    std::optional<RecordEntry> NextRecord();

    /// \brief Signal end-of-input. Throws if buffered bytes still contain a
    /// partial record or unparsed header.
    void Finish();

    [[nodiscard]] std::uint64_t RecordsEmitted() const noexcept { return recordsEmitted_; }

private:
    enum class State : std::uint8_t
    {
        NEED_MAGIC,
        NEED_L_TEXT,
        NEED_TEXT_BYTES,
        NEED_N_REF,
        NEED_REF_NAME_LEN,
        NEED_REF_NAME_AND_LEN,
        RECORDS,
    };

    struct BlockEntry
    {
        std::uint64_t startInStream;
        std::uint64_t fileOffset;
        std::uint32_t size;
    };

    std::vector<std::byte> buf_;
    std::uint64_t windowStart_{0};
    std::uint64_t cursor_{0};
    std::deque<BlockEntry> blocks_;
    State state_{State::NEED_MAGIC};
    std::uint32_t pendingTextLen_{0};
    std::int32_t pendingNRefRemaining_{0};
    std::uint32_t pendingRefNameLen_{0};
    std::uint64_t recordsEmitted_{0};

    [[nodiscard]] std::size_t BufferAvailable() const noexcept
    {
        return std::size(buf_) - (cursor_ - windowStart_);
    }

    [[nodiscard]] const std::byte* CursorPtr() const noexcept
    {
        return std::data(buf_) + (cursor_ - windowStart_);
    }

    void Advance(std::size_t n) { cursor_ += n; }

    void Compact()
    {
        const std::size_t consumed{cursor_ - windowStart_};
        if (consumed == 0U) {
            return;
        }
        buf_.erase(std::ranges::begin(buf_),
                   std::ranges::begin(buf_) + static_cast<std::ptrdiff_t>(consumed));
        windowStart_ = cursor_;
        while (!std::empty(blocks_) &&
               ((blocks_.front().startInStream + blocks_.front().size) <= windowStart_)) {
            blocks_.pop_front();
        }
    }

    [[nodiscard]] VirtualOffset VirtualOffsetForCursor() const
    {
        for (const BlockEntry& block : blocks_) {
            const std::uint64_t blockEnd{block.startInStream + block.size};
            if ((cursor_ >= block.startInStream) && (cursor_ < blockEnd)) {
                return VirtualOffset{block.fileOffset,
                                     static_cast<std::uint16_t>(cursor_ - block.startInStream)};
            }
        }
        throw std::runtime_error{"zmi-index: cursor not inside any tracked block"};
    }

    bool TryAdvanceHeader();
    std::optional<RecordEntry> TryEmitRecord();
};

inline void RecordParser::Feed(std::uint64_t fileOffset, std::span<const std::byte> data)
{
    Compact();
    const std::uint64_t blockStart{windowStart_ + std::size(buf_)};
    buf_.insert(std::ranges::end(buf_), std::ranges::begin(data), std::ranges::end(data));
    blocks_.push_back(
        BlockEntry{blockStart, fileOffset, static_cast<std::uint32_t>(std::size(data))});
}

inline bool RecordParser::TryAdvanceHeader()
{
    while (state_ != State::RECORDS) {
        switch (state_) {
            case State::NEED_MAGIC: {
                if (BufferAvailable() < 4U) {
                    return false;
                }
                const std::byte* p{CursorPtr()};
                if ((p[0] != std::byte{'B'}) || (p[1] != std::byte{'A'}) ||
                    (p[2] != std::byte{'M'}) || (p[3] != std::byte{0x01})) {
                    throw std::runtime_error{"zmi-index: BAM magic 'BAM\\1' not found"};
                }
                Advance(4U);
                state_ = State::NEED_L_TEXT;
                break;
            }
            case State::NEED_L_TEXT: {
                if (BufferAvailable() < 4U) {
                    return false;
                }
                pendingTextLen_ = ReadU32LE(CursorPtr());
                Advance(4U);
                state_ = State::NEED_TEXT_BYTES;
                break;
            }
            case State::NEED_TEXT_BYTES: {
                if (BufferAvailable() < pendingTextLen_) {
                    return false;
                }
                Advance(pendingTextLen_);
                state_ = State::NEED_N_REF;
                break;
            }
            case State::NEED_N_REF: {
                if (BufferAvailable() < 4U) {
                    return false;
                }
                pendingNRefRemaining_ = ReadI32LE(CursorPtr());
                Advance(4U);
                if (pendingNRefRemaining_ < 0) {
                    throw std::runtime_error{"zmi-index: BAM header n_ref negative"};
                }
                state_ = (pendingNRefRemaining_ == 0) ? State::RECORDS : State::NEED_REF_NAME_LEN;
                break;
            }
            case State::NEED_REF_NAME_LEN: {
                if (BufferAvailable() < 4U) {
                    return false;
                }
                pendingRefNameLen_ = ReadU32LE(CursorPtr());
                Advance(4U);
                state_ = State::NEED_REF_NAME_AND_LEN;
                break;
            }
            case State::NEED_REF_NAME_AND_LEN: {
                const std::size_t need{static_cast<std::size_t>(pendingRefNameLen_) + 4U};
                if (BufferAvailable() < need) {
                    return false;
                }
                Advance(need);
                --pendingNRefRemaining_;
                state_ = (pendingNRefRemaining_ == 0) ? State::RECORDS : State::NEED_REF_NAME_LEN;
                break;
            }
            case State::RECORDS:
                break;
        }
    }
    return true;
}

inline std::optional<RecordEntry> RecordParser::TryEmitRecord()
{
    if (BufferAvailable() < 4U) {
        return std::nullopt;
    }
    const std::uint32_t blockSize{ReadU32LE(CursorPtr())};
    if (BufferAvailable() < (4U + static_cast<std::size_t>(blockSize))) {
        return std::nullopt;
    }
    const VirtualOffset vo{VirtualOffsetForCursor()};
    const std::span<const std::byte> body{CursorPtr() + 4U, blockSize};
    const auto [rgId, zmw]{detail::ParseRecordIdentity(body)};
    Advance(4U + static_cast<std::size_t>(blockSize));
    ++recordsEmitted_;
    return RecordEntry{rgId, zmw, vo};
}

inline std::optional<RecordEntry> RecordParser::NextRecord()
{
    if (!TryAdvanceHeader()) {
        return std::nullopt;
    }
    return TryEmitRecord();
}

inline void RecordParser::Finish()
{
    if (state_ != State::RECORDS) {
        throw std::runtime_error{"zmi-index: input ended before BAM header completed"};
    }
    if (BufferAvailable() != 0U) {
        throw std::runtime_error{"zmi-index: input ended with partial record"};
    }
}

}  // namespace ZmiIndex
}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_TOOLS_ZMI_INDEX_RECORDPARSER_HPP
