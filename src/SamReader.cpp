#include <pbsamoa/io/SamReader.hpp>

#include <pbsamoa/core/CigarOp.hpp>
#include <pbsamoa/core/Tags.hpp>

#include <format>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <charconv>
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

std::int32_t ParseInt32(std::string_view text, std::string_view field)
{
    std::int32_t value{0};
    const auto [ptr, ec] =
        std::from_chars(std::data(text), std::data(text) + std::size(text), value);
    if (ec != std::errc{}) {
        throw std::runtime_error{
            std::format("SamReader: invalid integer '{}' for field {}", text, field)};
    }
    return value;
}

std::uint16_t ParseUInt16(std::string_view text, std::string_view field)
{
    std::uint16_t value{0};
    const auto [ptr, ec] =
        std::from_chars(std::data(text), std::data(text) + std::size(text), value);
    if (ec != std::errc{}) {
        throw std::runtime_error{
            std::format("SamReader: invalid integer '{}' for field {}", text, field)};
    }
    return value;
}

std::uint8_t ParseUInt8(std::string_view text, std::string_view field)
{
    // from_chars doesn't support uint8_t on all platforms, parse as uint16_t
    const std::uint16_t v{ParseUInt16(text, field)};
    return static_cast<std::uint8_t>(v);
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
            if (!headerResult.has_value()) {
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

std::optional<BamRecord> SamReader::ReadRecord()
{
    // Use buffered first alignment line if available
    if (impl_->hasBufferedLine) {
        impl_->hasBufferedLine = false;
        return ParseAlignmentLine(impl_->currentLine);
    }

    // Read next lines, skipping empty ones
    std::string line;
    while (std::getline(impl_->file, line)) {
        ++impl_->lineNumber;
        if (!std::empty(line)) {
            return ParseAlignmentLine(line);
        }
    }

    return std::nullopt;
}

BamRecord SamReader::ParseAlignmentLine(std::string_view line)
{
    const std::vector<std::string_view> fields{SplitOnTabs(line)};
    if (std::size(fields) < 11) {
        throw std::runtime_error{
            std::format("SamReader: line {}: expected at least 11 tab-separated fields, got {}",
                        impl_->lineNumber, std::size(fields))};
    }

    BamRecord record;

    // QNAME
    record.Name(std::string{fields[0]});

    // FLAG
    record.Flag(ParseUInt16(fields[1], "FLAG"));

    // RNAME → RefId
    if (fields[2] == "*") {
        record.RefId(-1);
    } else {
        record.RefId(impl_->header.ReferenceId(fields[2]));
    }

    // POS (1-based to 0-based)
    {
        const std::int32_t pos{ParseInt32(fields[3], "POS")};
        record.Pos((pos > 0) ? (pos - 1) : -1);
    }

    // MAPQ
    record.MapQ(ParseUInt8(fields[4], "MAPQ"));

    // CIGAR
    {
        auto cigar{ParseCigar(fields[5])};
        if (!cigar.has_value()) {
            throw std::runtime_error{
                std::format("SamReader: line {}: {}", impl_->lineNumber, cigar.error())};
        }
        record.Cigar(std::move(*cigar));
    }

    // RNEXT → NextRefId
    if (fields[6] == "*") {
        record.NextRefId(-1);
    } else if (fields[6] == "=") {
        record.NextRefId(record.RefId());
    } else {
        record.NextRefId(impl_->header.ReferenceId(fields[6]));
    }

    // PNEXT (1-based to 0-based)
    {
        const std::int32_t pnext{ParseInt32(fields[7], "PNEXT")};
        record.NextPos((pnext > 0) ? (pnext - 1) : -1);
    }

    // TLEN
    record.Tlen(ParseInt32(fields[8], "TLEN"));

    // SEQ
    if (fields[9] == "*") {
        record.Sequence(std::string{});
    } else {
        record.Sequence(std::string{fields[9]});
    }

    // QUAL
    if (fields[10] == "*") {
        record.Qualities(std::vector<std::uint8_t>{});
    } else {
        std::vector<std::uint8_t> quals;
        quals.reserve(std::size(fields[10]));
        for (const char c : fields[10]) {
            quals.push_back(c - 33);
        }
        record.Qualities(std::move(quals));
    }

    // Optional tags (fields 11+)
    TagMap tags;
    for (std::size_t i{11}; i < std::size(fields); ++i) {
        const auto parsed = ParseTagFromSam(fields[i]);
        if (parsed.has_value()) {
            tags.Set(parsed->first, parsed->second);
        }
    }
    if (tags.Size() > 0) {
        record.Tags(std::move(tags));
    }

    return record;
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
    if (!current_.has_value()) {
        reader_ = nullptr;
    }
}

const BamRecord& SamReader::RecordRange::Iterator::operator*() const { return *current_; }

const BamRecord* SamReader::RecordRange::Iterator::operator->() const { return &*current_; }

SamReader::RecordRange::Iterator& SamReader::RecordRange::Iterator::operator++()
{
    current_ = reader_->ReadRecord();
    if (!current_.has_value()) {
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
