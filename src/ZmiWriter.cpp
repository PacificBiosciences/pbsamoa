#include <pbsamoa/io/ZmiWriter.hpp>

#include "BinaryUtils.hpp"
#include "ZmiInternal.hpp"

#include <pbsamoa/core/Bgzf.hpp>

#include <array>

#include <cstdint>

namespace PacBio {
namespace Samoa {

struct ZmiWriter::Impl
{
    BgzfWriter bgzf;
    std::uint64_t numRecords{0};
    bool closed{false};

    explicit Impl(const std::filesystem::path& path, const ZmiWriterConfig& config)
        : bgzf{path, BgzfWriterConfig{
                         .CompressionLevel = 6,
                         .BgzfWorkers = 1,
                         .InputQueueCapacity = 256,
                         .BlocksPerBatch = 8,
                         .UseTempFile = config.UseTempFile,
                     }}
    {
        std::array<std::byte, detail::ZMI_HEADER_SIZE> header{};

        std::ranges::copy_n(std::begin(detail::ZMI_MAGIC), std::size(detail::ZMI_MAGIC),
                            std::begin(header));

        WriteLE(std::data(header) + 4, detail::ZMI_VERSION);

        WriteLE(std::data(header) + 8, detail::ZMI_ENTRY_SIZE_FIELD);

        bgzf.Write(header);
    }

    ~Impl() = default;
    Impl(const Impl&) = delete;
    Impl& operator=(const Impl&) = delete;
};

ZmiWriter::ZmiWriter(const std::filesystem::path& path, const ZmiWriterConfig& config)
    : impl_{std::make_unique<Impl>(path, config)}
{
}

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
}

}  // namespace Samoa
}  // namespace PacBio
