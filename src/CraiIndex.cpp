#include <pbsamoa/index/CraiIndex.hpp>

#include "BinaryUtils.hpp"
#include "CramInternal.hpp"

#include <libdeflate.h>

#include <algorithm>
#include <array>
#include <format>
#include <limits>
#include <memory>
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

std::runtime_error CraiGzipError(const std::filesystem::path& path)
{
    return std::runtime_error{std::format(
        "CraiIndex: failed to decompress gzip data in {} (truncated or invalid)", path.string())};
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
    std::size_t inputOffset{0};
    while (inputOffset < std::size(compressed)) {
        std::size_t outCapacity = std::max<std::size_t>(64, std::size(compressed) - inputOffset);
        bool memberDone{false};

        while (!memberDone) {
            memberOutput.resize(outCapacity);
            std::size_t consumed{0};
            std::size_t produced{0};
            const auto rc = libdeflate_gzip_decompress_ex(
                decompressor.get(), compressed.data() + inputOffset,
                std::size(compressed) - inputOffset, memberOutput.data(), std::size(memberOutput),
                &consumed, &produced);

            if (rc == LIBDEFLATE_SUCCESS) {
                if (consumed == 0) {
                    throw CraiGzipError(path);
                }
                decompressed.insert(std::ranges::end(decompressed),
                                    std::ranges::begin(memberOutput),
                                    std::ranges::begin(memberOutput) + produced);
                inputOffset += consumed;
                memberDone = true;
            } else if (rc == LIBDEFLATE_INSUFFICIENT_SPACE) {
                if (outCapacity >
                    std::numeric_limits<std::size_t>::max() / static_cast<std::size_t>(2)) {
                    throw CraiGzipError(path);
                }
                outCapacity *= 2;
            } else {
                throw CraiGzipError(path);
            }
        }
    }

    return decompressed;
}

std::array<std::string_view, CRAI_COLUMN_COUNT> ParseColumns(std::string_view line,
                                                             std::size_t lineNumber)
{
    std::array<std::string_view, CRAI_COLUMN_COUNT> columns{};
    std::size_t columnCount{0};
    std::size_t fieldStart{0};

    for (std::size_t i{0}; i <= std::size(line); ++i) {
        if (i == std::size(line) || line[i] == '\t') {
            if (columnCount < CRAI_COLUMN_COUNT) {
                columns[columnCount] = line.substr(fieldStart, i - fieldStart);
            }
            ++columnCount;
            fieldStart = i + 1;
        }
    }

    if (columnCount != CRAI_COLUMN_COUNT) {
        throw std::runtime_error{
            std::format("CraiIndex: line {} expected 6 columns, got {}", lineNumber, columnCount)};
    }

    return columns;
}

}  // namespace

CraiIndex CraiIndex::FromFile(const std::filesystem::path& path)
{
    const std::vector<std::byte> compressed = ReadAllBytes(path);
    const std::vector<std::byte> decompressed = DecompressGzipMembers(compressed, path);
    const std::string text{reinterpret_cast<const char*>(decompressed.data()),
                           std::size(decompressed)};

    CraiIndex index;
    std::size_t lineNumber{1};
    std::size_t lineStart{0};
    for (std::size_t i{0}; i <= std::size(text); ++i) {
        if (i == std::size(text) || text[i] == '\n') {
            std::string_view line{text.data() + lineStart, i - lineStart};
            if (!std::empty(line) && line.back() == '\r') {
                line.remove_suffix(1);
            }
            if (!std::empty(line)) {
                const auto columns = ParseColumns(line, lineNumber);
                const auto fieldCtx = [&](std::string_view field) {
                    return std::format("CraiIndex: line {} field {}", lineNumber, field);
                };
                index.AddEntry(CraiEntry{
                    .SequenceId = ParseInteger<std::int32_t>(columns[0], fieldCtx("sequence_id")),
                    .AlignmentStart =
                        ParseInteger<std::int64_t>(columns[1], fieldCtx("alignment_start")),
                    .AlignmentSpan =
                        ParseInteger<std::int64_t>(columns[2], fieldCtx("alignment_span")),
                    .ContainerOffset =
                        ParseInteger<std::int64_t>(columns[3], fieldCtx("container_offset")),
                    .SliceOffset = ParseInteger<std::int64_t>(columns[4], fieldCtx("slice_offset")),
                    .SliceSize = ParseInteger<std::int64_t>(columns[5], fieldCtx("slice_size")),
                });
            }
            lineStart = i + 1;
            ++lineNumber;
        }
    }

    return index;
}

const std::vector<CraiEntry>& CraiIndex::Entries() const { return entries_; }

std::vector<CraiEntry> CraiIndex::EntriesForReference(std::int32_t refId) const
{
    const auto it = referenceEntries_.find(refId);
    if (it == referenceEntries_.cend()) {
        return {};
    }

    std::vector<CraiEntry> result;
    result.reserve(std::size(it->second));
    for (const std::size_t index : it->second) {
        result.push_back(entries_.at(index));
    }
    return result;
}

void CraiIndex::AddEntry(CraiEntry entry)
{
    referenceEntries_[entry.SequenceId].push_back(std::size(entries_));
    entries_.push_back(std::move(entry));
}

}  // namespace Samoa
}  // namespace PacBio
