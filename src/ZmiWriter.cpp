#include <pbsamoa/io/ZmiWriter.hpp>

#include "ZmiInternal.hpp"

#include <pbsamoa/core/Bgzf.hpp>

#include <array>
#include <bit>
#include <type_traits>

#include <cstdint>

namespace PacBio {
namespace Samoa {

namespace {

template <typename T>
void WriteLE(std::byte* dst, T value)
{
    static_assert(std::is_integral_v<T>, "WriteLE requires an integral type");
    using UnsignedT = std::make_unsigned_t<T>;
    const UnsignedT bits{std::bit_cast<UnsignedT>(value)};
    for (std::size_t i{0}; i < sizeof(T); ++i) {
        dst[i] = static_cast<std::byte>((bits >> (i * 8U)) & static_cast<UnsignedT>(0xFFU));
    }
}

}  // namespace

struct ZmiWriter::Impl
{
    BgzfWriter bgzf;
    std::uint64_t numRecords{0};
    bool closed{false};

    explicit Impl(const std::filesystem::path& path)
        : bgzf{path, BgzfWriterConfig{
                         .CompressionLevel = 6,
                         .BgzfWorkers = 1,
                         .InputQueueCapacity = 256,
                         .BlocksPerBatch = 8,
                     }}
    {
        // Write the 64-byte header with numRecords = 0 as a placeholder.
        //
        // Design choice: numRecords is written as 0 (advisory). The reader
        // determines the actual record count by reading entries to EOF and
        // counting them. This avoids the complexity of patching the BGZF-
        // compressed header after close (which would require decompressing
        // the first block, patching, recompressing, and overwriting). The
        // numRecords field is reserved for future optimization where a reader
        // could pre-allocate if the value is nonzero.
        std::array<std::byte, detail::ZMI_HEADER_SIZE> header{};

        // Offset 0: magic "ZMI\1" (4 bytes)
        std::ranges::copy_n(std::begin(detail::ZMI_MAGIC), std::size(detail::ZMI_MAGIC),
                            std::begin(header));

        // Offset 4: version 0x010000 (4 bytes LE)
        WriteLE(std::data(header) + 4, detail::ZMI_VERSION);

        // Offset 8: entrySize = 16 (2 bytes LE)
        WriteLE(std::data(header) + 8, detail::ZMI_ENTRY_SIZE_FIELD);

        // Offset 10: flags = 0 (2 bytes LE, already zeroed)
        // Offset 12: numRecords = 0 (8 bytes LE, already zeroed)
        // Offset 20: reserved = zeros (44 bytes, already zeroed)

        bgzf.Write(header);
    }

    ~Impl() = default;
    Impl(const Impl&) = delete;
    Impl& operator=(const Impl&) = delete;
};

ZmiWriter::ZmiWriter(const std::filesystem::path& path) : impl_{std::make_unique<Impl>(path)} {}

ZmiWriter::~ZmiWriter()
{
    if ((impl_ != nullptr) && (!impl_->closed)) {
        Close();
    }
}

ZmiWriter::ZmiWriter(ZmiWriter&&) noexcept = default;
ZmiWriter& ZmiWriter::operator=(ZmiWriter&&) noexcept = default;

void ZmiWriter::AddRecord(std::int32_t rgId, std::int32_t zmw, std::int64_t virtualOffset)
{
    // Each entry is 16 bytes: rgId(4) + zmw(4) + virtualOffset(8), all
    // little-endian.
    std::array<std::byte, detail::ZMI_ENTRY_SIZE> entry{};
    WriteLE(std::data(entry), rgId);
    WriteLE(std::data(entry) + 4, zmw);
    WriteLE(std::data(entry) + 8, virtualOffset);

    impl_->bgzf.Write(entry);
    ++impl_->numRecords;
}

void ZmiWriter::Close()
{
    if (impl_->closed) {
        return;
    }

    impl_->bgzf.Close();
    impl_->closed = true;

    // Note: numRecords is tracked in impl_->numRecords for future use.
    // Currently the header contains numRecords=0 as advisory. The reader
    // counts entries by reading to EOF. See Impl constructor comment for
    // rationale on why we don't patch the header post-close.
}

}  // namespace Samoa
}  // namespace PacBio
