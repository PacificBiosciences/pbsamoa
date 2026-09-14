#include <pbsamoa/io/ZmiBamWriter.hpp>

#include "PathUtils.hpp"
#include "ZmwUtils.hpp"

#include <pbsamoa/core/RawRecord.hpp>
#include <pbsamoa/io/BamWriter.hpp>
#include <pbsamoa/io/ZmiWriter.hpp>

#include <cstddef>
#include <cstdint>
#include <exception>

namespace PacBio {
namespace Samoa {

struct ZmiBamWriter::Impl
{
    ZmiWriter zmi;
    BamWriter bam;
    bool closed{false};

    Impl(const std::filesystem::path& bamPath, const SamHeader& header, const BamWriterConfig& cfg)
        : zmi{SidecarPath(bamPath, ".zmi"), ZmiWriterConfig{.UseTempFile = cfg.UseTempFile}}
        , bam{bamPath, header, cfg,
              [this](std::int64_t virtualOffset, std::span<const std::byte> rawData) {
                  // This lambda reads the raw bytes directly. It does not build a
                  // RawRecord and a TagMap first, because it needs only RG, bc,
                  // and the read name. A full tag map would copy the PacBio ip/pw
                  // kinetics on the writer's IO thread, only to discard that data
                  // afterward.
                  BarcodePair barcodes;
                  const ZmwIdentity identity{detail::ParseRecordIdentity(rawData, &barcodes)};
                  zmi.AddRecord(identity.rgId, identity.zmw, virtualOffset, barcodes);
              }}
    {
    }

    Impl(const Impl&) = delete;
    Impl& operator=(const Impl&) = delete;
    ~Impl() = default;
};

ZmiBamWriter::ZmiBamWriter(const std::filesystem::path& bamPath, const SamHeader& header,
                           const BamWriterConfig& config)
    : impl_{std::make_unique<Impl>(bamPath, header, config)}
{
}

ZmiBamWriter::~ZmiBamWriter()
{
    if (impl_ && !impl_->closed) {
        Close();
    }
}

ZmiBamWriter::ZmiBamWriter(ZmiBamWriter&&) noexcept = default;
ZmiBamWriter& ZmiBamWriter::operator=(ZmiBamWriter&&) noexcept = default;

void ZmiBamWriter::Write(const BamRecord& record) { impl_->bam.Write(record); }

void ZmiBamWriter::Write(const RawRecordView& view) { impl_->bam.Write(view); }

void ZmiBamWriter::Write(const RawRecord& view) { impl_->bam.Write(view); }

void ZmiBamWriter::WriteBatch(const RawRecordBatch& batch) { impl_->bam.WriteBatch(batch); }

void ZmiBamWriter::Close()
{
    if (impl_->closed) {
        return;
    }
    impl_->closed = true;
    // BamWriter must close first (flushes remaining records + fires callbacks)
    std::exception_ptr error;
    try {
        impl_->bam.Close();
    } catch (...) {
        error = std::current_exception();
    }
    // Drain the index even when BAM serialization or an index callback failed.
    try {
        impl_->zmi.Close();
    } catch (...) {
        if (!error) {
            error = std::current_exception();
        }
    }
    if (error) {
        std::rethrow_exception(error);
    }
}

}  // namespace Samoa
}  // namespace PacBio
