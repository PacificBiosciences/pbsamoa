#ifndef PBSAMOA_CORE_SEQUENCE_HPP
#define PBSAMOA_CORE_SEQUENCE_HPP

#include <array>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {

/// \brief 256-entry lookup table: ASCII char -> 4-bit code.
/// Unmapped characters resolve to 15 (N).
inline constexpr std::array<std::uint8_t, 256> BASE_TO_CODE = [] {
    std::array<std::uint8_t, 256> t{};
    for (std::uint8_t& value : t) {
        value = 15;  // default: N
    }
    // =ACMGRSVTWYHKDBN -> 0-15
    constexpr std::string_view BASES{"=ACMGRSVTWYHKDBN"};
    for (std::uint8_t i{0}; i < 16; ++i) {
        t[static_cast<unsigned char>(BASES[i])] = i;
        // also handle lowercase
        if ((BASES[i] >= 'A') && (BASES[i] <= 'Z')) {
            t[static_cast<unsigned char>(BASES[i] + 32)] = i;
        }
    }
    return t;
}();

/// \brief 16-entry decode table: 4-bit code -> ASCII char.
inline constexpr std::array<char, 16> CODE_TO_BASE = {
    '=', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N',
};

/// \brief 256-entry lookup: packed byte -> two ASCII bases for fast decode.
inline constexpr std::array<std::array<char, 2>, 256> PACKED_TO_BASES = [] {
    std::array<std::array<char, 2>, 256> t{};
    for (std::size_t i{0}; i < std::size(t); ++i) {
        t[i][0] = CODE_TO_BASE[(i >> 4) & 0xF];
        t[i][1] = CODE_TO_BASE[i & 0xF];
    }
    return t;
}();

/// \brief Pack a text sequence into 4-bit encoding (2 bases per byte).
std::vector<std::byte> PackSequence(std::string_view seq);

/// \brief Pack a text sequence directly into a pre-allocated destination buffer.
///
/// Writes packed 4-bit bases into dest without allocating. Equivalent to
/// PackSequence but avoids the intermediate vector.
///
/// \param[in] seq      text sequence (ASCII bases)
/// \param[in] dest     destination buffer, must have at least (size(seq)+1)/2 bytes
void PackSequenceInto(std::string_view seq, std::byte* dest);

/// \brief Unpack a 4-bit encoded sequence to text.
std::string UnpackSequence(std::span<const std::byte> packed, std::uint32_t seqLength);

/// \brief Append decoded sequence directly to a string buffer.
void WriteSequenceTo(std::span<const std::byte> packed, std::uint32_t seqLength, std::string& out);

/// \brief Reverse complement a DNA sequence in-place.
void ReverseComplementInPlace(std::string& seq);

/// \brief Return the reverse complement of a DNA sequence.
std::string ReverseComplement(std::string_view seq);

/// \brief Non-owning view over 4-bit packed sequence bytes.
class SequenceView
{
public:
    SequenceView();
    SequenceView(std::span<const std::byte> packed, std::uint32_t length);

    char operator[](std::uint32_t i) const;
    std::uint32_t Size() const;
    std::string ToString() const;
    void WriteTo(std::string& out) const;

private:
    std::span<const std::byte> data_;
    std::uint32_t length_;
};

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_CORE_SEQUENCE_HPP
