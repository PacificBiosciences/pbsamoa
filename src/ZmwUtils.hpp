#ifndef PBSAMOA_ZMWUTILS_HPP
#define PBSAMOA_ZMWUTILS_HPP

#include <pbsamoa/core/Endian.hpp>
#include <pbsamoa/core/Tags.hpp>
#include <pbsamoa/index/ZmwIndex.hpp>

#include <algorithm>
#include <optional>
#include <span>
#include <stdexcept>
#include <string_view>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {

/// \brief Parse ZMW hole number from PacBio read name format "movie/<zmw>/...".
/// \returns hole number, or 0 on parse failure.
std::int32_t ParseZmwFromName(std::string_view name);

/// \brief Parse PacBio read-group ID string to its numeric int32 form.
/// \returns numeric ID, or 0 on parse failure.
std::int32_t ParseReadGroupId(std::string_view readGroupId);

/// \brief Parse the ZMW identity carried by a record name and tags.
/// \returns `(rgId, zmw)` with either field set to 0 when parsing fails.
ZmwIdentity ParseZmwIdentity(std::string_view name, const TagMap& tags);

namespace detail {

/// \brief Returns the length in bytes of the value part of a BAM aux field,
/// based on the type code.
///
/// \p data must point at the first byte of the value, after the 2-byte tag
/// key and the 1-byte type code. The function returns the consumed value
/// length, or std::nullopt if the type code is unknown or the buffer is
/// truncated.
[[nodiscard]] inline std::optional<std::size_t> AuxValueLength(std::byte typeByte,
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
            const auto it{std::ranges::find(data, std::byte{0})};
            if (it == std::ranges::end(data)) {
                return std::nullopt;
            }
            return static_cast<std::size_t>(it - std::ranges::begin(data)) + 1U;
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

/// \brief Searches a BAM aux block for the \p target tag and returns its
/// value as a string_view.
///
/// The result is meaningful only when \p target is a Z-typed tag, for
/// example RG. The function returns nullopt if it does not find the tag or
/// if the block is malformed.
[[nodiscard]] inline std::optional<std::string_view> FindZTag(std::span<const std::byte> aux,
                                                              TagKey target)
{
    std::size_t i{0};
    while ((i + 3U) <= std::size(aux)) {
        const TagKey key{static_cast<char>(aux[i]), static_cast<char>(aux[i + 1U])};
        const std::byte typeByte{aux[i + 2U]};
        const std::span<const std::byte> rest{aux.subspan(i + 3U)};
        const std::optional<std::size_t> valueLen{AuxValueLength(typeByte, rest)};
        // AuxValueLength derives fixed-width and B-array lengths from the encoding
        // alone. A truncated value therefore yields a length past the end of
        // `rest`. The code rejects that value here so that the loop does not
        // advance `i` past the buffer.
        if (!valueLen || (*valueLen > std::size(rest))) {
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

/// \brief Parses the ZMW identity of a BAM record body.
///
/// \p recordBody holds the bytes that follow the 4-byte block_size prefix.
/// Its size equals the block_size value.
///
/// The function reads the name and the RG tag directly. It skips every
/// other tag by length and does not decode it. This technique means the
/// function never decodes the PacBio ip/pw kinetics only to obtain a read
/// group value.
///
/// \throws std::runtime_error on a structurally malformed record.
inline ZmwIdentity ParseRecordIdentity(std::span<const std::byte> recordBody)
{
    constexpr std::size_t FIXED_FIELDS_SIZE{32U};
    if (std::size(recordBody) < FIXED_FIELDS_SIZE) {
        throw std::runtime_error{"ParseRecordIdentity: BAM record smaller than fixed header"};
    }

    const std::uint8_t nameLen{static_cast<std::uint8_t>(recordBody[8])};
    const std::uint16_t cigarOpCount{ReadU16LE(std::data(recordBody) + 12U)};
    // The code declares seqLength as std::size_t, even though ReadU32LE reads
    // the underlying value as a 32-bit unsigned integer. This choice matters.
    // If l_seq equals UINT32_MAX, a 32-bit computation of (seqLength + 1U)
    // wraps to 0. Then auxOffset points back inside the record instead of past
    // it. This behavior needs std::size_t to stay 64 bits wide. Endian.hpp
    // enforces that width with a static_assert.
    const std::size_t seqLength{ReadU32LE(std::data(recordBody) + 16U)};

    if (nameLen == 0U) {
        throw std::runtime_error{"ParseRecordIdentity: BAM record has zero l_read_name"};
    }

    const std::size_t nameOffset{FIXED_FIELDS_SIZE};
    const std::size_t cigarOffset{nameOffset + nameLen};
    const std::size_t seqOffset{cigarOffset + (std::size_t{4U} * cigarOpCount)};
    const std::size_t qualOffset{seqOffset + ((seqLength + 1U) / 2U)};
    const std::size_t auxOffset{qualOffset + seqLength};

    if (auxOffset > std::size(recordBody)) {
        throw std::runtime_error{
            "ParseRecordIdentity: BAM record variable-length fields exceed block_size"};
    }
    if (recordBody[nameOffset + nameLen - 1U] != std::byte{0}) {
        throw std::runtime_error{"ParseRecordIdentity: BAM record name not NUL-terminated"};
    }

    const std::string_view name{reinterpret_cast<const char*>(std::data(recordBody) + nameOffset),
                                static_cast<std::size_t>(nameLen) - 1U};
    const std::span<const std::byte> aux{recordBody.subspan(auxOffset)};

    const std::optional<std::string_view> rgText{FindZTag(aux, RG_TAG)};
    return ZmwIdentity{rgText ? ParseReadGroupId(*rgText) : 0, ParseZmwFromName(name)};
}

}  // namespace detail

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_ZMWUTILS_HPP
