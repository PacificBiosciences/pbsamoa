#include <pbsamoa/io/SamReader.hpp>

#include "BinaryUtils.hpp"

#include <pbsamoa/core/CigarOp.hpp>
#include <pbsamoa/core/Tags.hpp>

#include <format>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <cstdint>

namespace PacBio {
namespace Samoa {

namespace {

/// \brief Split a string_view on tab characters into a vector of fields.
std::vector<std::string_view> SplitOnTabs(std::string_view line)
{
    std::vector<std::string_view> fields;
    while (true) {
        const std::size_t pos{line.find('\t')};
        if (pos == std::string_view::npos) {
            fields.push_back(line);
            break;
        }
        fields.push_back(line.substr(0, pos));
        line.remove_prefix(pos + 1);
    }
    return fields;
}

/// \brief Parse a uint8 from text via uint16 (from_chars may not support
/// uint8_t).
std::uint8_t ParseUInt8(std::string_view text, std::string_view context)
{
    const std::uint16_t v{ParseInteger<std::uint16_t>(text, context)};
    if (v > std::numeric_limits<std::uint8_t>::max()) {
        throw std::runtime_error{
            std::format("{}: value out of range (0-255): '{}'", context, text)};
    }
    return static_cast<std::uint8_t>(v);
}

std::int32_t ToZeroBasedSamPosition(const std::int32_t oneBasedPos)
{
    if (oneBasedPos <= 0) {
        return -1;
    }
    return oneBasedPos - 1;
}

std::int32_t ParseReferenceId(std::string_view refName, const SamHeader& header,
                              std::uint64_t lineNumber, std::string_view fieldName)
{
    if (refName == "*") {
        return -1;
    }

    const std::int32_t refId{header.ReferenceId(refName)};
    if ((header.NumReferences() > 0) && (refId < 0)) {
        throw std::runtime_error{
            std::format("SamReader: line {}: unknown {} '{}'", lineNumber, fieldName, refName)};
    }
    return refId;
}

}  // namespace

struct SamReader::Impl
{
    std::ifstream file;
    SamHeader header;
    std::string currentLine;
    std::uint64_t lineNumber{0};
    bool hasBufferedLine{false};

    explicit Impl(const std::filesystem::path& path)
    {
        file.open(path);
        if (!file.is_open()) {
            throw std::runtime_error{std::format("SamReader: cannot open file: {}", path.string())};
        }

        // Read header lines (starting with @)
        std::string headerText;
        std::string line;
        while (std::getline(file, line)) {
            ++lineNumber;
            if ((!std::empty(line)) && (line[0] == '@')) {
                headerText += line;
                headerText += '\n';
            } else {
                // First non-header line — buffer it
                currentLine = std::move(line);
                hasBufferedLine = !std::empty(currentLine);
                break;
            }
        }

        if (!std::empty(headerText)) {
            auto headerResult{SamHeader::FromText(headerText)};
            if (!headerResult) {
                throw std::runtime_error{std::move(headerResult.error())};
            }
            header = std::move(*headerResult);
        }
    }

    ~Impl() = default;
    Impl(const Impl&) = delete;
    Impl& operator=(const Impl&) = delete;
};

SamReader::SamReader(const std::filesystem::path& path) : impl_{std::make_unique<Impl>(path)} {}

SamReader::~SamReader() = default;
SamReader::SamReader(SamReader&&) noexcept = default;
SamReader& SamReader::operator=(SamReader&&) noexcept = default;

const SamHeader& SamReader::Header() const { return impl_->header; }

namespace {

BamRecord ParseAlignmentLine(std::string_view line, const SamHeader& header,
                             std::uint64_t lineNumber)
{
    const std::vector<std::string_view> fields{SplitOnTabs(line)};
    if (std::size(fields) < 11) {
        throw std::runtime_error{
            std::format("SamReader: line {}: expected at least 11 tab-separated fields, got {}",
                        lineNumber, std::size(fields))};
    }

    BamRecord record;

    // QNAME
    record.Name(std::string{fields[0]});

    // FLAG
    record.Flag(ParseInteger<std::uint16_t>(fields[1], "FLAG"));

    // RNAME → RefId
    record.RefId(ParseReferenceId(fields[2], header, lineNumber, "RNAME"));

    // POS (1-based to 0-based)
    const std::int32_t pos{ParseInteger<std::int32_t>(fields[3], "POS")};
    record.Pos(ToZeroBasedSamPosition(pos));

    // MAPQ
    record.MapQ(ParseUInt8(fields[4], "MAPQ"));

    // CIGAR
    auto cigar{ParseCigar(fields[5])};
    if (!cigar) {
        throw std::runtime_error{std::format("SamReader: line {}: {}", lineNumber, cigar.error())};
    }
    record.Cigar(std::move(*cigar));

    // RNEXT → NextRefId
    if (fields[6] == "=") {
        record.NextRefId(record.RefId());
    } else {
        record.NextRefId(ParseReferenceId(fields[6], header, lineNumber, "RNEXT"));
    }

    // PNEXT (1-based to 0-based)
    const std::int32_t pnext{ParseInteger<std::int32_t>(fields[7], "PNEXT")};
    record.NextPos(ToZeroBasedSamPosition(pnext));

    // TLEN
    record.Tlen(ParseInteger<std::int32_t>(fields[8], "TLEN"));

    // SEQ
    const bool seqMissing{fields[9] == "*"};
    record.Sequence(seqMissing ? std::string{} : std::string{fields[9]});

    // QUAL
    const bool qualMissing{fields[10] == "*"};
    if (qualMissing) {
        record.Qualities(std::vector<std::uint8_t>{});
    } else {
        if (seqMissing) {
            throw std::runtime_error{
                std::format("SamReader: line {}: QUAL provided while SEQ is '*'", lineNumber)};
        }

        std::vector<std::uint8_t> quals;
        quals.reserve(std::size(fields[10]));
        for (const char c : fields[10]) {
            if ((c < 33) || (c > 126)) {
                throw std::runtime_error{
                    std::format("SamReader: line {}: invalid QUAL character", lineNumber)};
            }
            quals.push_back(c - 33);
        }
        if (std::size(quals) != std::size(fields[9])) {
            throw std::runtime_error{
                std::format("SamReader: line {}: SEQ and QUAL lengths differ", lineNumber)};
        }
        record.Qualities(std::move(quals));
    }

    // Optional tags (fields 11+)
    TagMap tags;
    for (std::size_t i{11}; i < std::size(fields); ++i) {
        const auto parsedTag{ParseTagFromSam(fields[i])};
        if (!parsedTag) {
            throw std::runtime_error{
                std::format("SamReader: line {}: invalid tag '{}'", lineNumber, fields[i])};
        }
        tags.Set(parsedTag->first, parsedTag->second);
    }
    if (tags.Size() > 0) {
        record.Tags(std::move(tags));
    }

    return record;
}

}  // namespace

std::optional<BamRecord> SamReader::ReadRecord()
{
    // Use buffered first alignment line if available
    if (impl_->hasBufferedLine) {
        impl_->hasBufferedLine = false;
        return ParseAlignmentLine(impl_->currentLine, impl_->header, impl_->lineNumber);
    }

    // Read next lines, skipping empty ones
    std::string line;
    while (std::getline(impl_->file, line)) {
        ++impl_->lineNumber;
        if (!std::empty(line)) {
            return ParseAlignmentLine(line, impl_->header, impl_->lineNumber);
        }
    }

    return std::nullopt;
}

// --- RecordRange ---

SamReader::RecordRange::RecordRange(SamReader* reader) : reader_{reader} {}

SamReader::RecordRange::Iterator SamReader::RecordRange::begin() { return Iterator{reader_}; }

SamReader::RecordRange::Iterator SamReader::RecordRange::end() { return Iterator{}; }

// --- Iterator ---

SamReader::RecordRange::Iterator::Iterator() = default;

SamReader::RecordRange::Iterator::Iterator(SamReader* reader) : reader_{reader}
{
    current_ = reader_->ReadRecord();
    if (!current_) {
        reader_ = nullptr;
    }
}

const BamRecord& SamReader::RecordRange::Iterator::operator*() const { return *current_; }

const BamRecord* SamReader::RecordRange::Iterator::operator->() const { return &*current_; }

SamReader::RecordRange::Iterator& SamReader::RecordRange::Iterator::operator++()
{
    current_ = reader_->ReadRecord();
    if (!current_) {
        reader_ = nullptr;
    }
    return *this;
}

void SamReader::RecordRange::Iterator::operator++(int) { ++(*this); }

bool SamReader::RecordRange::Iterator::operator==(const Iterator& other) const
{
    return reader_ == other.reader_;
}

SamReader::RecordRange SamReader::Records() { return RecordRange{this}; }

}  // namespace Samoa
}  // namespace PacBio
