#include <pbsamoa/core/CigarOp.hpp>

#include <array>
#include <charconv>
#include <expected>
#include <format>

namespace {

void AppendBinRange(std::vector<std::uint16_t>& bins, std::int32_t begin, std::int32_t end,
                    std::int32_t offset)
{
    for (std::int32_t bin{offset + begin}; bin <= (offset + end); ++bin) {
        bins.push_back(static_cast<std::uint16_t>(bin));
    }
}

}  // namespace

namespace PacBio {
namespace Samoa {

std::expected<std::vector<CigarOp>, std::string> ParseCigar(std::string_view text)
{
    if (std::empty(text) || text == "*") {
        return std::vector<CigarOp>{};
    }

    std::vector<CigarOp> result;
    const char* pos{std::data(text)};
    const char* const end{std::data(text) + std::size(text)};

    while (pos < end) {
        std::uint32_t length{0};
        const std::from_chars_result parseResult{std::from_chars(pos, end, length)};
        if (parseResult.ec != std::errc{}) {
            return std::unexpected{"Invalid CIGAR: bad integer"};
        }
        if (parseResult.ptr == end) {
            return std::unexpected{"Invalid CIGAR: missing op character"};
        }
        const char* const opPos{parseResult.ptr};

        const std::optional<CigarOpType> opType{CharToCigarOp(*opPos)};
        if (!opType.has_value()) {
            return std::unexpected{std::format("Invalid CIGAR op: {}", *opPos)};
        }

        if (length == 0) {
            return std::unexpected{"Invalid CIGAR: operation length must be >= 1"};
        }
        result.emplace_back(*opType, length);
        pos = opPos + 1;
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
        const std::to_chars_result toCharsResult{
            std::to_chars(std::data(numBuf), std::data(numBuf) + std::size(numBuf), op.Length())};
        out.append(std::data(numBuf), toCharsResult.ptr);
        out += CigarOpToChar(op.Type());
    }
}

std::vector<std::uint16_t> Reg2Bins(std::int32_t beg, std::int32_t end)
{
    const std::int32_t inclusiveEnd{end - 1};
    std::vector<std::uint16_t> bins;
    bins.reserve(32);
    bins.push_back(0);
    AppendBinRange(bins, beg >> 26, inclusiveEnd >> 26, 1);
    AppendBinRange(bins, beg >> 23, inclusiveEnd >> 23, 9);
    AppendBinRange(bins, beg >> 20, inclusiveEnd >> 20, 73);
    AppendBinRange(bins, beg >> 17, inclusiveEnd >> 17, 585);
    AppendBinRange(bins, beg >> 14, inclusiveEnd >> 14, 4681);
    return bins;
}

}  // namespace Samoa
}  // namespace PacBio
