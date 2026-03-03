#include <pbsamoa/core/Sequence.hpp>

#include <array>
#include <stdexcept>

namespace PacBio {
namespace Samoa {

namespace {
inline constexpr std::array<char, 256> COMPLEMENT_TABLE = [] {
    std::array<char, 256> t{};
    for (int i{0}; i < 256; ++i) {
        t[i] = 'N';
    }
    t['A'] = 'T';
    t['T'] = 'A';
    t['C'] = 'G';
    t['G'] = 'C';
    t['N'] = 'N';
    t['a'] = 't';
    t['t'] = 'a';
    t['c'] = 'g';
    t['g'] = 'c';
    t['n'] = 'n';
    // IUPAC ambiguity codes
    t['M'] = 'K';
    t['K'] = 'M';
    t['R'] = 'Y';
    t['Y'] = 'R';
    t['S'] = 'S';
    t['W'] = 'W';
    t['B'] = 'V';
    t['V'] = 'B';
    t['D'] = 'H';
    t['H'] = 'D';
    return t;
}();

void DecodePackedInto(std::span<const std::byte> packed, std::uint32_t seqLength, char* dest)
{
    const std::uint32_t fullBytes{seqLength / 2};
    for (std::uint32_t i{0}; i < fullBytes; ++i) {
        const auto& pair{PACKED_TO_BASES[static_cast<std::uint8_t>(packed[i])]};
        dest[2 * i] = pair[0];
        dest[2 * i + 1] = pair[1];
    }
    if ((seqLength % 2) != 0) {
        dest[seqLength - 1] =
            CODE_TO_BASE[(static_cast<std::uint8_t>(packed[fullBytes]) >> 4) & 0xF];
    }
}
}  // namespace

std::vector<std::byte> PackSequence(std::string_view seq)
{
    const std::size_t packedSize{(std::size(seq) + 1) / 2};
    std::vector<std::byte> result(packedSize);

    for (std::size_t i{0}; i < std::size(seq); ++i) {
        const std::uint8_t code{BASE_TO_CODE[static_cast<unsigned char>(seq[i])]};
        if ((i % 2) == 0) {
            result[i / 2] = static_cast<std::byte>(code << 4);
        } else {
            result[i / 2] |= static_cast<std::byte>(code);
        }
    }
    return result;
}

void PackSequenceInto(std::string_view seq, std::byte* dest)
{
    const std::size_t len{std::size(seq)};
    const std::size_t fullPairs{len / 2};
    for (std::size_t i{0}; i < fullPairs; ++i) {
        const std::uint8_t hi{BASE_TO_CODE[static_cast<unsigned char>(seq[2 * i])]};
        const std::uint8_t lo{BASE_TO_CODE[static_cast<unsigned char>(seq[2 * i + 1])]};
        dest[i] = static_cast<std::byte>((hi << 4) | lo);
    }
    if ((len % 2) != 0) {
        dest[fullPairs] =
            static_cast<std::byte>(BASE_TO_CODE[static_cast<unsigned char>(seq[len - 1])] << 4);
    }
}

std::string UnpackSequence(std::span<const std::byte> packed, std::uint32_t seqLength)
{
    std::string result(seqLength, '\0');
    char* dest{std::data(result)};

    const std::uint32_t fullBytes{seqLength / 2};
    for (std::uint32_t i{0}; i < fullBytes; ++i) {
        const auto& pair{PACKED_TO_BASES[static_cast<std::uint8_t>(packed[i])]};
        dest[2 * i] = pair[0];
        dest[2 * i + 1] = pair[1];
    }
    if ((seqLength % 2) != 0) {
        dest[seqLength - 1] =
            CODE_TO_BASE[(static_cast<std::uint8_t>(packed[fullBytes]) >> 4) & 0xF];
    }
    return result;
}

void WriteSequenceTo(std::span<const std::byte> packed, std::uint32_t seqLength, std::string& out)
{
    const std::size_t startPos{std::size(out)};
    out.resize(startPos + seqLength);
    DecodePackedInto(packed, seqLength, std::data(out) + startPos);
}

SequenceView::SequenceView() : data_{}, length_{0} {}

SequenceView::SequenceView(std::span<const std::byte> packed, std::uint32_t length)
    : data_{packed}, length_{length}
{
}

char SequenceView::operator[](std::uint32_t i) const
{
    if (i >= length_) {
        throw std::out_of_range{"SequenceView::operator[]: index out of range"};
    }
    const std::uint8_t byte{static_cast<std::uint8_t>(data_[i / 2])};
    const std::uint8_t code{((i % 2) == 0) ? static_cast<std::uint8_t>(byte >> 4)
                                           : static_cast<std::uint8_t>(byte & 0x0FU)};
    return CODE_TO_BASE[code];
}

std::uint32_t SequenceView::Size() const { return length_; }

std::string SequenceView::ToString() const { return UnpackSequence(data_, length_); }

void SequenceView::WriteTo(std::string& out) const { WriteSequenceTo(data_, length_, out); }

void ReverseComplementInPlace(std::string& seq)
{
    const std::size_t n{std::size(seq)};
    if (n == 0) {
        return;
    }
    for (std::size_t i{0}; i < (n + 1) / 2; ++i) {
        const std::size_t j{n - 1 - i};
        const char ci{COMPLEMENT_TABLE[static_cast<unsigned char>(seq[i])]};
        const char cj{COMPLEMENT_TABLE[static_cast<unsigned char>(seq[j])]};
        seq[i] = cj;
        seq[j] = ci;
    }
}

std::string ReverseComplement(std::string_view seq)
{
    std::string result{seq};
    ReverseComplementInPlace(result);
    return result;
}

}  // namespace Samoa
}  // namespace PacBio
