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
    table['='] = '=';  // '=' (match-to-reference) is self-complementary; keep it lossless
    // IUPAC ambiguity codes (upper + lower case, matching htslib's
    // case-preserving revcomp)
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
    table['m'] = 'k';
    table['k'] = 'm';
    table['r'] = 'y';
    table['y'] = 'r';
    table['s'] = 's';
    table['w'] = 'w';
    table['b'] = 'v';
    table['v'] = 'b';
    table['d'] = 'h';
    table['h'] = 'd';
    return table;
}

constexpr std::array<char, 256> COMPLEMENT_TABLE = MakeComplementTable();

// The packed byte of two bases, by (first base << 8) | second base: one lookup for each byte,
// not two. PackSequenceInto takes a large part of the CRAM decode of HiFi reads.
std::array<std::uint8_t, 65536> MakePairToPackedTable()
{
    std::array<std::uint8_t, 65536> table{};
    for (std::size_t first = 0; first < 256; ++first) {
        for (std::size_t second = 0; second < 256; ++second) {
            table[(first << 8) | second] = (BASE_TO_CODE[first] << 4) | BASE_TO_CODE[second];
        }
    }
    return table;
}

// Built at run time, once: as a constexpr table, its 65536 checked subscripts go over the
// constexpr step limit of Clang when _GLIBCXX_ASSERTIONS is on.
const std::array<std::uint8_t, 65536>& PairToPacked()
{
    static const std::array<std::uint8_t, 65536> table{MakePairToPackedTable()};
    return table;
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

void PackSequenceInto(std::string_view seq, std::byte* dest)
{
    const std::size_t len{std::size(seq)};
    std::size_t packedIndex{0};
    std::size_t i{0};
    const std::array<std::uint8_t, 65536>& pairToPacked{PairToPacked()};
    for (; (i + 1) < len; i += 2, ++packedIndex) {
        const std::size_t pair{(std::size_t{static_cast<unsigned char>(seq[i])} << 8) |
                               static_cast<unsigned char>(seq[i + 1])};
        dest[packedIndex] = std::byte{pairToPacked[pair]};
    }
    if (i < len) {
        dest[packedIndex] =
            static_cast<std::byte>(BASE_TO_CODE[static_cast<unsigned char>(seq[i])] << 4);
    }
}

std::vector<std::byte> PackSequence(std::string_view seq)
{
    const std::size_t packedSize{(std::size(seq) + 1) / 2};
    std::vector<std::byte> result(packedSize);
    PackSequenceInto(seq, std::data(result));
    return result;
}

std::string UnpackSequence(std::span<const std::byte> packed, std::uint32_t seqLength)
{
    std::string result;
    result.resize_and_overwrite(seqLength, [&](char* buf, std::size_t size) {
        DecodePackedInto(packed, seqLength, buf);
        return size;
    });
    return result;
}

void WriteSequenceTo(std::span<const std::byte> packed, std::uint32_t seqLength, std::string& out)
{
    const std::size_t startPos{std::size(out)};
    out.resize_and_overwrite(startPos + seqLength, [&](char* buf, std::size_t size) {
        DecodePackedInto(packed, seqLength, buf + startPos);
        return size;
    });
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
