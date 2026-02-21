#include <pbsamoa/core/CigarOp.hpp>

#include <array>
#include <charconv>
#include <format>
#include <stdexcept>

namespace PacBio {
namespace Samoa {

std::vector<CigarOp> ParseCigar(std::string_view text)
{
    if (std::empty(text) || text == "*") {
        return {};
    }

    std::vector<CigarOp> result;
    const char* pos{std::data(text)};
    const char* end{std::data(text) + std::size(text)};

    while (pos < end) {
        std::uint32_t length{0};
        const std::from_chars_result parseResult{std::from_chars(pos, end, length)};
        const char* ptr{parseResult.ptr};
        const std::errc ec{parseResult.ec};
        if ((ec != std::errc{}) || (ptr == end)) {
            if (ec != std::errc{}) {
                throw std::runtime_error{"Invalid CIGAR: bad integer"};
            }
            throw std::runtime_error{"Invalid CIGAR: missing op character"};
        }

        const std::optional<CigarOpType> opType{CharToCigarOp(*ptr)};
        if (!opType.has_value()) {
            throw std::runtime_error{std::string{"Invalid CIGAR op: "} + *ptr};
        }

        if (length == 0) {
            throw std::runtime_error{"Invalid CIGAR: operation length must be >= 1"};
        }
        result.emplace_back(*opType, length);
        pos = ptr + 1;
    }

    return result;
}

std::string CigarToString(CigarView cigar)
{
    if (std::empty(cigar)) {
        return "*";
    }

    std::string result;
    for (const CigarOp& op : cigar) {
        std::format_to(std::back_inserter(result), "{}", op.Length());
        result += CigarOpToChar(op.Type());
    }
    return result;
}

void WriteCigarTo(CigarView cigar, std::string& out)
{
    if (std::empty(cigar)) {
        out += '*';
        return;
    }

    std::array<char, 16> numBuf{};
    for (const CigarOp& op : cigar) {
        const auto [ptr, ec]{
            std::to_chars(std::data(numBuf), std::data(numBuf) + std::size(numBuf), op.Length())};
        out.append(std::data(numBuf), ptr);
        out += CigarOpToChar(op.Type());
    }
}

std::vector<std::uint16_t> Reg2Bins(std::int32_t beg, std::int32_t end)
{
    std::vector<std::uint16_t> bins;
    bins.reserve(32);
    --end;
    bins.push_back(0);
    for (std::int32_t k = 1 + (beg >> 26); k <= 1 + (end >> 26); ++k) {
        bins.push_back(k);
    }
    for (std::int32_t k = 9 + (beg >> 23); k <= 9 + (end >> 23); ++k) {
        bins.push_back(k);
    }
    for (std::int32_t k = 73 + (beg >> 20); k <= 73 + (end >> 20); ++k) {
        bins.push_back(k);
    }
    for (std::int32_t k = 585 + (beg >> 17); k <= 585 + (end >> 17); ++k) {
        bins.push_back(k);
    }
    for (std::int32_t k = 4681 + (beg >> 14); k <= 4681 + (end >> 14); ++k) {
        bins.push_back(k);
    }
    return bins;
}

}  // namespace Samoa
}  // namespace PacBio
