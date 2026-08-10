#ifndef PBSAMOA_TOOLS_ZMI_INDEX_RECORDPARSER_HPP
#define PBSAMOA_TOOLS_ZMI_INDEX_RECORDPARSER_HPP

#include <pbsamoa/core/Bgzf.hpp>
#include <pbsamoa/core/Endian.hpp>
#include <pbsamoa/core/Tags.hpp>

#include "../../ZmwUtils.hpp"

#include <deque>
#include <optional>
#include <span>
#include <stdexcept>
#include <string_view>
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

// AuxValueLength, FindZTag, and ParseRecordIdentity now live in ZmwUtils.hpp. The ZMI
// writer needs the same walker functions. An unqualified detail:: reference here
// resolves to the enclosing PacBio::Samoa::detail namespace.

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
