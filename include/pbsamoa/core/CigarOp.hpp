#ifndef PBSAMOA_CORE_CIGAROP_HPP
#define PBSAMOA_CORE_CIGAROP_HPP

#include <array>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <cstdint>

namespace PacBio {
namespace Samoa {

/// \brief CIGAR operation type codes matching BAM encoding (spec §4.2).
enum class CigarOpType : std::uint8_t
{
    M = 0,   // alignment match (sequence match or mismatch)
    I = 1,   // insertion to the reference
    D = 2,   // deletion from the reference
    N = 3,   // skipped region from the reference
    S = 4,   // soft clipping
    H = 5,   // hard clipping
    P = 6,   // padding
    EQ = 7,  // sequence match
    X = 8,   // sequence mismatch
};

/// \brief Single CIGAR operation: 32-bit value matching BAM binary layout.
/// Stored as `op_len << 4 | op` so spans of CigarOp can alias BAM data.
class CigarOp
{
public:
    constexpr CigarOp();
    constexpr explicit CigarOp(std::uint32_t rawValue);
    constexpr CigarOp(CigarOpType type, std::uint32_t length);

    constexpr CigarOpType Type() const;
    constexpr std::uint32_t Length() const;
    constexpr std::uint32_t RawValue() const;

    constexpr bool operator==(const CigarOp&) const = default;

private:
    std::uint32_t data_;
};

/// \brief Non-owning view over a CIGAR operation sequence.
using CigarView = std::span<const CigarOp>;

// --- free functions ---

constexpr char CigarOpToChar(CigarOpType op);
constexpr std::optional<CigarOpType> CharToCigarOp(char c);
constexpr bool ConsumesQuery(CigarOpType op);
constexpr bool ConsumesReference(CigarOpType op);

/// \brief Reference-consuming length of a CIGAR string.
constexpr std::int64_t ReferenceLength(CigarView cigar);

/// \brief Query-consuming length of a CIGAR string.
constexpr std::int64_t QueryLength(CigarView cigar);

/// \brief Parse CIGAR from SAM text (e.g., "8M2I4M1D3M"). "*" returns empty.
std::vector<CigarOp> ParseCigar(std::string_view text);

/// \brief Serialize CIGAR operations to SAM text. Empty returns "*".
std::string CigarToString(CigarView cigar);

/// \brief Append CIGAR text directly to a string buffer. Empty appends "*".
void WriteCigarTo(CigarView cigar, std::string& out);

/// \brief Compute BAI index bin for a 0-based half-open interval [beg, end).
constexpr std::uint16_t Reg2Bin(std::int32_t beg, std::int32_t end);

/// \brief Compute all BAI bins overlapping a 0-based half-open interval [beg, end).
std::vector<std::uint16_t> Reg2Bins(std::int32_t beg, std::int32_t end);

// --- inline constexpr definitions ---

constexpr CigarOp::CigarOp() : data_{0} {}

constexpr CigarOp::CigarOp(std::uint32_t rawValue) : data_{rawValue} {}

constexpr CigarOp::CigarOp(CigarOpType type, std::uint32_t length)
    : data_{(length << 4) | static_cast<std::uint32_t>(type)}
{
}

constexpr CigarOpType CigarOp::Type() const { return static_cast<CigarOpType>(data_ & 0xFU); }

constexpr std::uint32_t CigarOp::Length() const { return data_ >> 4; }

constexpr std::uint32_t CigarOp::RawValue() const { return data_; }

constexpr char CigarOpToChar(CigarOpType op)
{
    constexpr std::array<char, 9> TABLE{'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X'};
    const std::uint8_t idx{std::to_underlying(op)};
    return (idx < std::size(TABLE)) ? TABLE[idx] : '?';
}

constexpr std::optional<CigarOpType> CharToCigarOp(char c)
{
    switch (c) {
        case 'M':
            return CigarOpType::M;
        case 'I':
            return CigarOpType::I;
        case 'D':
            return CigarOpType::D;
        case 'N':
            return CigarOpType::N;
        case 'S':
            return CigarOpType::S;
        case 'H':
            return CigarOpType::H;
        case 'P':
            return CigarOpType::P;
        case '=':
            return CigarOpType::EQ;
        case 'X':
            return CigarOpType::X;
        default:
            return std::nullopt;
    }
}

constexpr bool ConsumesQuery(CigarOpType op)
{
    // M, I, S, =, X consume query
    switch (op) {
        case CigarOpType::M:
        case CigarOpType::I:
        case CigarOpType::S:
        case CigarOpType::EQ:
        case CigarOpType::X:
            return true;
        default:
            return false;
    }
}

constexpr bool ConsumesReference(CigarOpType op)
{
    // M, D, N, =, X consume reference
    switch (op) {
        case CigarOpType::M:
        case CigarOpType::D:
        case CigarOpType::N:
        case CigarOpType::EQ:
        case CigarOpType::X:
            return true;
        default:
            return false;
    }
}

constexpr std::int64_t ReferenceLength(CigarView cigar)
{
    std::int64_t len{0};
    for (const CigarOp& op : cigar) {
        if (ConsumesReference(op.Type())) {
            len += op.Length();
        }
    }
    return len;
}

constexpr std::int64_t QueryLength(CigarView cigar)
{
    std::int64_t len{0};
    for (const CigarOp& op : cigar) {
        if (ConsumesQuery(op.Type())) {
            len += op.Length();
        }
    }
    return len;
}

constexpr std::uint16_t Reg2Bin(std::int32_t beg, std::int32_t end)
{
    --end;
    if ((beg >> 14) == (end >> 14)) {
        return ((1 << 15) - 1) / 7 + (beg >> 14);
    }
    if ((beg >> 17) == (end >> 17)) {
        return ((1 << 12) - 1) / 7 + (beg >> 17);
    }
    if ((beg >> 20) == (end >> 20)) {
        return ((1 << 9) - 1) / 7 + (beg >> 20);
    }
    if ((beg >> 23) == (end >> 23)) {
        return ((1 << 6) - 1) / 7 + (beg >> 23);
    }
    if ((beg >> 26) == (end >> 26)) {
        return ((1 << 3) - 1) / 7 + (beg >> 26);
    }
    return 0;
}

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_CORE_CIGAROP_HPP
