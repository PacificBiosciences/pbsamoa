#include <pbsamoa/io/FastaReference.hpp>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>
#include <charconv>
#include <cstring>
#include <format>
#include <fstream>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <unordered_map>
#include <utility>
#include <vector>

#include <cstdint>

namespace PacBio {
namespace Samoa {

namespace {

// One sequence of the file. A regular sequence is described like a .fai line; an irregular one
// (lines of uneven length, blank lines, spaces) is held as text.
struct FastaSequence
{
    std::int64_t Length{0};
    std::int64_t Offset{0};     // of the first base in the file
    std::int64_t LineBases{0};  // bases per full line
    std::int64_t LineBytes{0};  // bytes per full line, with the line end
    std::optional<std::string> Text;
};

}  // namespace

struct FastaFile::Data
{
    std::filesystem::path Path;
    const char* Map{nullptr};
    std::size_t Size{0};
    std::vector<std::string> Names;  // in file order
    std::vector<FastaSequence> Sequences;
    std::unordered_map<std::string, std::size_t> Index;  // name -> first sequence of that name

    Data() = default;
    Data(const Data&) = delete;
    Data& operator=(const Data&) = delete;

    ~Data()
    {
        if (Map != nullptr) {
            ::munmap(const_cast<char*>(Map), Size);
        }
    }

    void Add(std::string name, FastaSequence sequence)
    {
        Index.try_emplace(name, std::size(Sequences));
        Names.push_back(std::move(name));
        Sequences.push_back(std::move(sequence));
    }

    // The bases [start, start + length) of sequence i (inside it); a view of the file when they
    // lie on one line, else copied to buffer without the line ends.
    std::string_view Bases(const std::size_t i, const std::int64_t start, const std::int64_t length,
                           std::string& buffer) const
    {
        const FastaSequence& s{Sequences[i]};
        if (s.Text) {
            return std::string_view{*s.Text}.substr(static_cast<std::size_t>(start),
                                                    static_cast<std::size_t>(length));
        }
        if (length <= 0) {
            return {};
        }
        std::int64_t line{start / s.LineBases};
        std::int64_t column{start % s.LineBases};
        if (column + length <= s.LineBases) {
            return {Map + s.Offset + (line * s.LineBytes) + column,
                    static_cast<std::size_t>(length)};
        }
        buffer.resize(static_cast<std::size_t>(length));
        std::int64_t copied{0};
        while (copied < length) {
            const std::int64_t n{std::min(s.LineBases - column, length - copied)};
            std::memcpy(std::data(buffer) + copied, Map + s.Offset + (line * s.LineBytes) + column,
                        static_cast<std::size_t>(n));
            copied += n;
            ++line;
            column = 0;
        }
        return buffer;
    }
};

namespace {

bool IsGzip(const char* map, const std::size_t size)
{
    return (size >= 2) && (static_cast<unsigned char>(map[0]) == 0x1f) &&
           (static_cast<unsigned char>(map[1]) == 0x8b);
}

// A view of the lines of the mapped file: each without its '\n', and with the offset of the
// next line.
class LineScanner
{
public:
    LineScanner(const char* map, const std::size_t size) : map_{map}, size_{size} {}

    bool Next(std::string_view& line, std::size_t& lineStart)
    {
        if (pos_ >= size_) {
            return false;
        }
        lineStart = pos_;
        const void* nl{std::memchr(map_ + pos_, '\n', size_ - pos_)};
        const std::size_t end{(nl == nullptr)
                                  ? size_
                                  : static_cast<std::size_t>(static_cast<const char*>(nl) - map_)};
        line = std::string_view{map_ + pos_, end - pos_};
        pos_ = (nl == nullptr) ? size_ : end + 1;
        return true;
    }

    std::size_t Position() const { return pos_; }

private:
    const char* map_;
    std::size_t size_;
    std::size_t pos_{0};
};

std::string HeaderName(std::string_view line)
{
    line.remove_prefix(1);
    return std::string{line.substr(0, line.find_first_of(" \t\r"))};
}

// The bases of the lines of one sequence, without whitespace; for irregular sequences.
std::string CollectText(const char* map, const std::size_t begin, const std::size_t end)
{
    std::string text;
    text.reserve(end - begin);
    for (std::size_t i = begin; i < end; ++i) {
        const char c{map[i]};
        if ((c != '\n') && (c != '\r') && (c != ' ') && (c != '\t') && (c != '\v') && (c != '\f')) {
            text.push_back(c);
        }
    }
    return text;
}

// Indexes the mapped file in one pass over its lines.
void ScanFasta(FastaFile::Data& data)
{
    LineScanner lines{data.Map, data.Size};
    std::string_view line;
    std::size_t lineStart{0};

    std::optional<std::string> name;
    FastaSequence seq;
    std::size_t seqBegin{0};
    std::int64_t lineCount{0};
    bool lastLineShort{false};  // a line shorter than LineBases was seen: it must be the last
    bool blankSeen{false};
    bool irregular{false};

    const auto finish = [&](const std::size_t seqEnd) {
        if (!name) {
            return;
        }
        if (irregular) {
            seq.Text = CollectText(data.Map, seqBegin, seqEnd);
            seq.Length = std::ssize(*seq.Text);
        }
        data.Add(std::move(*name), std::move(seq));
        name.reset();
    };

    while (lines.Next(line, lineStart)) {
        if (!std::empty(line) && (line.front() == '>')) {
            finish(lineStart);
            name = HeaderName(line);
            seq = FastaSequence{};
            seqBegin = lines.Position();
            lineCount = 0;
            lastLineShort = false;
            blankSeen = false;
            irregular = false;
            continue;
        }
        if (!name || irregular) {
            continue;
        }
        const bool crlf{!std::empty(line) && (line.back() == '\r')};
        const std::string_view bases{crlf ? line.substr(0, std::size(line) - 1) : line};
        if (std::empty(bases)) {
            blankSeen = true;
            continue;
        }
        const auto lineBytes{static_cast<std::int64_t>(lines.Position() - lineStart)};
        const bool endsInNewline{lines.Position() > lineStart + std::size(line)};
        const auto n{std::ssize(bases)};
        // The last line of the file may lack its '\n'; every other full line needs LineBytes.
        if (blankSeen || lastLineShort ||
            (bases.find_first_of(" \t\v\f") != std::string_view::npos) ||
            ((lineCount > 0) && ((n > seq.LineBases) || ((n == seq.LineBases) && endsInNewline &&
                                                         (lineBytes != seq.LineBytes))))) {
            irregular = true;
            continue;
        }
        if (lineCount == 0) {
            seq.Offset = static_cast<std::int64_t>(lineStart);
            seq.LineBases = n;
            seq.LineBytes = lineBytes;
        } else if (n < seq.LineBases) {
            lastLineShort = true;
        }
        seq.Length += n;
        ++lineCount;
    }
    finish(data.Size);
}

std::int64_t ParseFaiField(const std::string_view field, const std::filesystem::path& fai)
{
    std::int64_t value{0};
    const auto [ptr, ec] =
        std::from_chars(std::data(field), std::data(field) + std::size(field), value);
    if ((ec != std::errc{}) || (ptr != std::data(field) + std::size(field)) || (value < 0)) {
        throw std::runtime_error{
            std::format("FASTA index {}: bad field '{}'", fai.string(), field)};
    }
    return value;
}

// Reads <fasta>.fai: name, length, offset, line bases, line bytes (and more columns for FASTQ).
void ReadFai(FastaFile::Data& data, const std::filesystem::path& fai)
{
    std::ifstream in{fai};
    if (!in) {
        throw std::runtime_error{std::format("cannot open FASTA index {}", fai.string())};
    }
    std::string text;
    while (std::getline(in, text)) {
        if (std::empty(text)) {
            continue;
        }
        std::vector<std::string_view> fields;
        std::string_view rest{text};
        while (true) {
            const auto tab{rest.find('\t')};
            fields.push_back(rest.substr(0, tab));
            if (tab == std::string_view::npos) {
                break;
            }
            rest.remove_prefix(tab + 1);
        }
        if (std::size(fields) < 5) {
            throw std::runtime_error{
                std::format("FASTA index {}: expected 5 columns in '{}'", fai.string(), text)};
        }
        FastaSequence seq{
            .Length = ParseFaiField(fields[1], fai),
            .Offset = ParseFaiField(fields[2], fai),
            .LineBases = ParseFaiField(fields[3], fai),
            .LineBytes = ParseFaiField(fields[4], fai),
        };
        if (seq.Length > 0) {
            const std::int64_t last{seq.Length - 1};
            const std::int64_t end{
                seq.Offset + ((last / std::max<std::int64_t>(seq.LineBases, 1)) * seq.LineBytes) +
                (last % std::max<std::int64_t>(seq.LineBases, 1)) + 1};
            if ((seq.LineBases == 0) || (seq.LineBytes < seq.LineBases) ||
                (end > static_cast<std::int64_t>(data.Size))) {
                throw std::runtime_error{std::format(
                    "FASTA index {} does not fit {} (sequence {}); rebuild it with samtools faidx",
                    fai.string(), data.Path.string(), fields[0])};
            }
        }
        data.Add(std::string{fields[0]}, std::move(seq));
    }
}

// The reference of the @SQ lines of a header: an index into the sequences of the file.
class FastaFileReference final : public CramReference
{
public:
    FastaFileReference(std::shared_ptr<const FastaFile::Data> data,
                       std::vector<std::optional<std::size_t>> sequences)
        : data_{std::move(data)}, sequences_{std::move(sequences)}
    {
    }

    std::int32_t NumSequences() const override
    {
        return static_cast<std::int32_t>(std::size(sequences_));
    }

    std::int64_t Length(const std::int32_t refId) const override
    {
        const auto& i{sequences_[static_cast<std::size_t>(refId)]};
        return i ? data_->Sequences[*i].Length : 0;
    }

    std::string_view Bases(const std::int32_t refId, const std::int64_t start,
                           const std::int64_t length, std::string& buffer) const override
    {
        const auto& i{sequences_[static_cast<std::size_t>(refId)]};
        if (!i) {
            return {};
        }
        return data_->Bases(*i, start, length, buffer);
    }

private:
    std::shared_ptr<const FastaFile::Data> data_;
    std::vector<std::optional<std::size_t>> sequences_;
};

}  // namespace

FastaFile::FastaFile(const std::filesystem::path& path)
{
    auto data{std::make_shared<Data>()};
    data->Path = path;

    const int fd{::open(path.c_str(), O_RDONLY | O_CLOEXEC)};
    if (fd < 0) {
        throw std::runtime_error{std::format("cannot open reference FASTA {}", path.string())};
    }
    struct stat st{};
    if (::fstat(fd, &st) != 0) {
        ::close(fd);
        throw std::runtime_error{std::format("cannot stat reference FASTA {}", path.string())};
    }
    data->Size = static_cast<std::size_t>(st.st_size);
    if (data->Size > 0) {
        void* map{::mmap(nullptr, data->Size, PROT_READ, MAP_PRIVATE, fd, 0)};
        if (map == MAP_FAILED) {
            ::close(fd);
            throw std::runtime_error{std::format("cannot map reference FASTA {}", path.string())};
        }
        data->Map = static_cast<const char*>(map);
    }
    ::close(fd);

    if ((data->Map != nullptr) && !IsGzip(data->Map, data->Size)) {
        const std::filesystem::path fai{path.string() + ".fai"};
        std::error_code ec;
        if (std::filesystem::is_regular_file(fai, ec)) {
            ReadFai(*data, fai);
        } else {
            ScanFasta(*data);
        }
    }
    data_ = std::move(data);
}

std::size_t FastaFile::NumSequences() const { return std::size(data_->Sequences); }

std::shared_ptr<const CramReference> FastaFile::Reference(const SamHeader& header,
                                                          const MissingFastaSequence missing) const
{
    std::vector<std::optional<std::size_t>> sequences;
    sequences.reserve(std::size(header.ReferenceSequences()));
    for (const ReferenceSequence& sq : header.ReferenceSequences()) {
        const auto it{data_->Index.find(std::string{sq.Name()})};
        if (it == std::end(data_->Index)) {
            if (missing == MissingFastaSequence::THROW) {
                throw std::runtime_error{
                    std::format("reference sequence {} is not in the FASTA file", sq.Name())};
            }
            sequences.emplace_back();
        } else {
            sequences.emplace_back(it->second);
        }
    }
    return std::make_shared<const FastaFileReference>(data_, std::move(sequences));
}

}  // namespace Samoa
}  // namespace PacBio
