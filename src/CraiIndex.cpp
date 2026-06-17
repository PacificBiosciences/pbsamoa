#include <pbsamoa/index/CraiIndex.hpp>

#include "BinaryUtils.hpp"
#include "LibdeflateUtils.hpp"

#include <libdeflate.h>

#include <algorithm>
#include <array>
#include <format>
#include <limits>
#include <ranges>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace PacBio {
namespace Samoa {

namespace {

constexpr std::size_t CRAI_COLUMN_COUNT{6};

struct GzipMemberResult
{
    std::size_t consumed;
    std::size_t produced;
};

std::runtime_error CraiGzipError(const std::filesystem::path& path)
{
    return std::runtime_error{std::format(
        "CraiIndex: failed to decompress gzip data in {} (truncated or invalid)", path.string())};
}

GzipMemberResult DecompressGzipMember(const LibdeflateDecompressorPtr& decompressor,
                                      std::span<const std::byte> compressed,
                                      std::size_t inputOffset, std::vector<std::byte>& memberOutput,
                                      const std::filesystem::path& path)
{
    std::size_t outCapacity{std::max<std::size_t>(64, std::size(compressed) - inputOffset)};

    while (true) {
        memberOutput.resize(outCapacity);

        std::size_t consumed{0};
        std::size_t produced{0};

        const enum libdeflate_result rc {
            libdeflate_gzip_decompress_ex(decompressor.get(), compressed.data() + inputOffset,
                                          std::size(compressed) - inputOffset, memberOutput.data(),
                                          std::size(memberOutput), &consumed, &produced)
        };

        if (rc == LIBDEFLATE_SUCCESS) {
            if (consumed == 0) {
                throw CraiGzipError(path);
            }
            return GzipMemberResult{consumed, produced};
        }
        if (rc == LIBDEFLATE_INSUFFICIENT_SPACE) {
            if (outCapacity >
                std::numeric_limits<std::size_t>::max() / static_cast<std::size_t>(2)) {
                throw CraiGzipError(path);
            }
            outCapacity *= 2;
            continue;
        }
        throw CraiGzipError(path);
    }
}

std::vector<std::byte> DecompressGzipMembers(std::span<const std::byte> compressed,
                                             const std::filesystem::path& path)
{
    const LibdeflateDecompressorPtr decompressor{libdeflate_alloc_decompressor()};
    if (!decompressor) {
        throw CraiGzipError(path);
    }

    std::vector<std::byte> decompressed;
    std::vector<std::byte> memberOutput;
    const std::size_t compressedSize{std::size(compressed)};
    std::size_t inputOffset{0};
    while (inputOffset < compressedSize) {
        const GzipMemberResult member{
            DecompressGzipMember(decompressor, compressed, inputOffset, memberOutput, path)};
        decompressed.insert(std::ranges::end(decompressed), std::ranges::begin(memberOutput),
                            std::ranges::begin(memberOutput) + member.produced);
        inputOffset += member.consumed;
    }

    return decompressed;
}

std::string FieldContext(std::size_t lineNumber, std::string_view field)
{
    return std::format("CraiIndex: line {} field {}", lineNumber, field);
}

std::array<std::string_view, CRAI_COLUMN_COUNT> ParseColumns(std::string_view line,
                                                             std::size_t lineNumber)
{
    std::array<std::string_view, CRAI_COLUMN_COUNT> columns{};
    std::size_t columnCount{0};
    std::size_t fieldStart{0};
    const std::size_t lineSize{std::size(line)};

    for (std::size_t i{0}; i <= lineSize; ++i) {
        if (i == lineSize || line[i] == '\t') {
            if (columnCount < CRAI_COLUMN_COUNT) {
                columns[columnCount] = line.substr(fieldStart, i - fieldStart);
            }
            ++columnCount;
            fieldStart = i + 1;
        }
    }

    if (columnCount != CRAI_COLUMN_COUNT) {
        throw std::runtime_error{std::format("CraiIndex: line {} expected {} columns, got {}",
                                             lineNumber, CRAI_COLUMN_COUNT, columnCount)};
    }

    return columns;
}

CraiEntry ParseEntry(std::span<const std::string_view, CRAI_COLUMN_COUNT> columns,
                     std::size_t lineNumber)
{
    return CraiEntry{
        .SequenceId =
            ParseInteger<std::int32_t>(columns[0], FieldContext(lineNumber, "sequence_id")),
        .AlignmentStart =
            ParseInteger<std::int64_t>(columns[1], FieldContext(lineNumber, "alignment_start")),
        .AlignmentSpan =
            ParseInteger<std::int64_t>(columns[2], FieldContext(lineNumber, "alignment_span")),
        .ContainerOffset =
            ParseInteger<std::int64_t>(columns[3], FieldContext(lineNumber, "container_offset")),
        .SliceOffset =
            ParseInteger<std::int64_t>(columns[4], FieldContext(lineNumber, "slice_offset")),
        .SliceSize = ParseInteger<std::int64_t>(columns[5], FieldContext(lineNumber, "slice_size")),
    };
}

}  // namespace

CraiIndex CraiIndex::FromFile(const std::filesystem::path& path)
{
    const std::vector<std::byte> compressed = ReadAllBytes(path);
    const std::vector<std::byte> decompressed = DecompressGzipMembers(compressed, path);
    const std::string_view text{reinterpret_cast<const char*>(decompressed.data()),
                                std::size(decompressed)};
    const std::size_t textSize{std::size(text)};

    CraiIndex index;
    std::size_t lineNumber{1};
    std::size_t lineStart{0};
    for (std::size_t i{0}; i <= textSize; ++i) {
        if (i == textSize || text[i] == '\n') {
            std::string_view line{text.data() + lineStart, i - lineStart};
            if (!std::empty(line) && line.back() == '\r') {
                line.remove_suffix(1);
            }
            if (!std::empty(line)) {
                const auto columns{ParseColumns(line, lineNumber)};
                index.AddEntry(ParseEntry(columns, lineNumber));
            }
            lineStart = i + 1;
            ++lineNumber;
        }
    }

    return index;
}

std::vector<CraiEntry> CraiIndex::EntriesForReference(std::int32_t refId) const
{
    const auto it{entries_.find(refId)};
    if (it == std::end(entries_)) {
        return {};
    }
    return it->second;
}

void CraiIndex::AddEntry(CraiEntry entry)
{
    entries_[entry.SequenceId].push_back(std::move(entry));
}

}  // namespace Samoa
}  // namespace PacBio
