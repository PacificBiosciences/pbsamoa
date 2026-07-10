#include <pbsamoa/core/CigarOp.hpp>

#include <array>
#include <charconv>
#include <expected>
#include <format>

namespace {

struct BinLevel
{
    std::int32_t shift;
    std::int32_t offset;
};

constexpr std::array<BinLevel, 5> BIN_LEVELS{{
    {26, 1},
    {23, 9},
    {20, 73},
    {17, 585},
    {14, 4681},
}};

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
        const auto [ptr, ec]{std::from_chars(pos, end, length)};
        if (ec != std::errc{}) {
            return std::unexpected{"Invalid CIGAR: bad integer"};
        }
        if (ptr == end) {
            return std::unexpected{"Invalid CIGAR: missing op character"};
        }
        const char* const opPos{ptr};

        const std::optional<CigarOpType> opType{CharToCigarOp(*opPos)};
        if (!opType) {
            return std::unexpected{std::format("Invalid CIGAR op: {}", *opPos)};
        }

        // The BAM CIGAR field packs the length into 28 bits (len << 4 | op); a larger value
        // would corrupt the op code. std::from_chars already rejects > UINT32_MAX. A zero
        // length is grammar-valid and accepted by htslib, so it is not rejected here.
        if (length > 0x0FFFFFFFu) {
            return std::unexpected{"Invalid CIGAR: operation length out of range [0, 268435455]"};
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
    // The binning scheme spans at most 2^29 bp; clamp the query end like htslib's reg2bins
    // (hts.c) so an out-of-range end cannot generate bins outside the table.
    constexpr std::int32_t MAX_BIN_COORD{1 << 29};
    if (end > MAX_BIN_COORD) {
        end = MAX_BIN_COORD;
    }
    if (beg >= end) {
        return {};
    }
    const std::int32_t inclusiveEnd{end - 1};
    std::vector<std::uint16_t> bins;
    bins.reserve(32);
    bins.push_back(0);
    for (const BinLevel& level : BIN_LEVELS) {
        const std::int32_t beginBin{beg >> level.shift};
        const std::int32_t endBin{inclusiveEnd >> level.shift};
        for (std::int32_t bin{level.offset + beginBin}; bin <= (level.offset + endBin); ++bin) {
            bins.push_back(static_cast<std::uint16_t>(bin));
        }
    }
    return bins;
}

CigarStats ComputeCigarStats(CigarView cigar)
{
    // Raw counts, matching pbcopper Data::CigarOpsCalculator. MatchBases counts
    // both = and M; X is mismatch; indels accumulate both bases and event counts.
    std::int32_t matchBases{0};
    std::int32_t mismatchBases{0};
    std::int32_t insertionBases{0};
    std::int32_t deletionBases{0};
    std::int32_t insertionEvents{0};
    std::int32_t deletionEvents{0};

    for (const CigarOp& op : cigar) {
        const std::int32_t len{static_cast<std::int32_t>(op.Length())};
        switch (op.Type()) {
            case CigarOpType::EQ:
            case CigarOpType::M:
                matchBases += len;
                break;
            case CigarOpType::X:
                mismatchBases += len;
                break;
            case CigarOpType::I:
                insertionBases += len;
                ++insertionEvents;
                break;
            case CigarOpType::D:
                deletionBases += len;
                ++deletionEvents;
                break;
            case CigarOpType::N:
            case CigarOpType::S:
            case CigarOpType::H:
            case CigarOpType::P:
                break;
        }
    }

    CigarStats stats;
    stats.Mismatches = mismatchBases;
    stats.Insertions = insertionBases;
    stats.Deletions = deletionBases;
    stats.InsertionEvents = insertionEvents;
    stats.DeletionEvents = deletionEvents;
    stats.NumAlignedBases = matchBases + insertionBases + mismatchBases;
    stats.Span = stats.NumAlignedBases;

    const std::int32_t denomIdentity{matchBases + mismatchBases + deletionBases + insertionBases};
    if (denomIdentity > 0) {
        stats.Identity = (100.0 * matchBases) / denomIdentity;
    }

    const std::int32_t denomGapComp{matchBases + mismatchBases + deletionEvents + insertionEvents};
    if (denomGapComp > 0) {
        stats.IdentityGapComp = (100.0 * matchBases) / denomGapComp;
    }

    return stats;
}

}  // namespace Samoa
}  // namespace PacBio
