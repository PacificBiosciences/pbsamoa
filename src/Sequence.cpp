#include <pbsamoa/core/Sequence.hpp>

#include <array>
#include <stdexcept>

namespace PacBio {
namespace Samoa {

namespace {
constexpr std::array<char, 256> MakeComplementTable()
{
    std::array<char, 256> table{};
    table.fill('N');
    table['A'] = 'T';
    table['T'] = 'A';
    table['C'] = 'G';
    table['G'] = 'C';
    table['N'] = 'N';
    table['a'] = 't';
    table['t'] = 'a';
    table['c'] = 'g';
    table['g'] = 'c';
    table['n'] = 'n';
    // IUPAC ambiguity codes
    table['M'] = 'K';
    table['K'] = 'M';
    table['R'] = 'Y';
    table['Y'] = 'R';
    table['S'] = 'S';
    table['W'] = 'W';
    table['B'] = 'V';
    table['V'] = 'B';
    table['D'] = 'H';
    table['H'] = 'D';
    return table;
}

inline constexpr std::array<char, 256> COMPLEMENT_TABLE = MakeComplementTable();

void PackSequenceIntoImpl(std::string_view seq, std::byte* dest)
{
    const std::size_t len{std::size(seq)};
    std::size_t packedIndex{0};
    std::size_t i{0};
    for (; (i + 1) < len; i += 2, ++packedIndex) {
        const std::uint8_t hi{BASE_TO_CODE[static_cast<unsigned char>(seq[i])]};
        const std::uint8_t lo{BASE_TO_CODE[static_cast<unsigned char>(seq[i + 1])]};
        dest[packedIndex] = static_cast<std::byte>((hi << 4) | lo);
    }
    if (i < len) {
        dest[packedIndex] =
            static_cast<std::byte>(BASE_TO_CODE[static_cast<unsigned char>(seq[i])] << 4);
    }
}

void DecodePackedInto(std::span<const std::byte> packed, std::uint32_t seqLength, char* dest)
{
    const std::uint32_t fullBytes{seqLength / 2};
    for (std::uint32_t i{0}; i < fullBytes; ++i) {
        const std::array<char, 2>& decodedBases{
            PACKED_TO_BASES[static_cast<std::uint8_t>(packed[i])]};
        dest[2 * i] = decodedBases[0];
        dest[2 * i + 1] = decodedBases[1];
    }
    if ((seqLength % 2) != 0) {
        const std::uint8_t packedByte{static_cast<std::uint8_t>(packed[fullBytes])};
        dest[seqLength - 1] = CODE_TO_BASE[packedByte >> 4];
    }
}
}  // namespace

std::vector<std::byte> PackSequence(std::string_view seq)
{
    const std::size_t len{std::size(seq)};
    const std::size_t packedSize{(len + 1) / 2};
    std::vector<std::byte> result(packedSize);
    PackSequenceIntoImpl(seq, std::data(result));
    return result;
}

void PackSequenceInto(std::string_view seq, std::byte* dest) { PackSequenceIntoImpl(seq, dest); }

std::string UnpackSequence(std::span<const std::byte> packed, std::uint32_t seqLength)
{
    std::string result(seqLength, '\0');
    DecodePackedInto(packed, seqLength, std::data(result));
    return result;
}

void WriteSequenceTo(std::span<const std::byte> packed, std::uint32_t seqLength, std::string& out)
{
    const std::size_t startPos{std::size(out)};
    out.resize(startPos + seqLength);
    DecodePackedInto(packed, seqLength, std::data(out) + startPos);
}

SequenceView::SequenceView() = default;

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
    if ((i % 2) == 0) {
        return CODE_TO_BASE[static_cast<std::uint8_t>(byte >> 4)];
    }
    return CODE_TO_BASE[static_cast<std::uint8_t>(byte & 0x0FU)];
}

std::uint32_t SequenceView::Size() const { return length_; }

std::string SequenceView::ToString() const { return UnpackSequence(data_, length_); }

void SequenceView::WriteTo(std::string& out) const { WriteSequenceTo(data_, length_, out); }

void ReverseComplementInPlace(std::string& seq)
{
    const std::size_t n{std::size(seq)};
    const std::size_t half{(n + 1) / 2};
    for (std::size_t i{0}; i < half; ++i) {
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
