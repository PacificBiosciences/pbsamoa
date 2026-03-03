#include <pbsamoa/cram/CramCompression.hpp>

#include "CramInternal.hpp"

#include <libdeflate.h>
#ifdef PBSAMOA_HAVE_BZIP2
#include <bzlib.h>
#endif
#include <htscodecs/arith_dynamic.h>
#include <htscodecs/fqzcomp_qual.h>
#include <htscodecs/rANS_static.h>
#include <htscodecs/rANS_static4x16.h>
#include <htscodecs/tokenise_name3.h>
#ifdef PBSAMOA_HAVE_LZMA
#include <lzma.h>
#endif

#include <algorithm>
#include <format>
#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <utility>

#include <cstdlib>

namespace PacBio {
namespace Samoa {

std::vector<std::byte> CramFqzcompCompress(std::span<const std::byte> data,
                                           std::span<const std::uint32_t> recordLengths,
                                           std::span<const std::uint32_t> recordFlags);

namespace {

std::vector<std::byte> CopyAndFree(unsigned char* ptr, std::size_t size, const char* context)
{
    if (!ptr) {
        throw std::runtime_error(std::format("{}: codec returned null output", context));
    }
    std::vector<std::byte> output(size);
    if (size > 0) {
        std::ranges::copy_n(reinterpret_cast<const std::byte*>(ptr), size, output.data());
    }
    std::free(ptr);
    return output;
}

std::vector<std::byte> CopyAndFree(char* ptr, std::size_t size, const char* context)
{
    return CopyAndFree(reinterpret_cast<unsigned char*>(ptr), size, context);
}

unsigned char* NonConstBytes(std::span<const std::byte> data)
{
    return reinterpret_cast<unsigned char*>(const_cast<std::byte*>(data.data()));
}

std::vector<std::byte> CramRans4x8Decompress(std::span<const std::byte> data, std::size_t rawSize)
{
    if (rawSize == 0) {
        return {};
    }
    if (std::size(data) > std::numeric_limits<unsigned int>::max() ||
        rawSize > std::numeric_limits<unsigned int>::max()) {
        throw std::runtime_error("CramRans4x8Decompress: input/output size too large");
    }

    unsigned int outSize = 0;
    auto* out =
        rans_uncompress(NonConstBytes(data), static_cast<unsigned int>(std::size(data)), &outSize);
    if (outSize != rawSize) {
        std::free(out);
        throw std::runtime_error("CramRans4x8Decompress: decompressed size mismatch");
    }
    return CopyAndFree(out, outSize, "CramRans4x8Decompress");
}

std::vector<std::byte> CramRans4x8Compress(std::span<const std::byte> data)
{
    if (std::empty(data)) {
        return {};
    }
    if (std::size(data) > std::numeric_limits<unsigned int>::max()) {
        throw std::runtime_error("CramRans4x8Compress: input size too large");
    }

    unsigned int outSize = 0;
    auto* out =
        rans_compress(NonConstBytes(data), static_cast<unsigned int>(std::size(data)), &outSize, 1);
    return CopyAndFree(out, outSize, "CramRans4x8Compress");
}

std::vector<std::byte> CramRans4x16Decompress(std::span<const std::byte> data, std::size_t rawSize)
{
    if (rawSize == 0) {
        return {};
    }
    if (std::size(data) > std::numeric_limits<unsigned int>::max() ||
        rawSize > std::numeric_limits<unsigned int>::max()) {
        throw std::runtime_error("CramRans4x16Decompress: input/output size too large");
    }

    unsigned int outSize = 0;
    auto* out = rans_uncompress_4x16(NonConstBytes(data),
                                     static_cast<unsigned int>(std::size(data)), &outSize);
    if (outSize != rawSize) {
        std::free(out);
        throw std::runtime_error("CramRans4x16Decompress: decompressed size mismatch");
    }
    return CopyAndFree(out, outSize, "CramRans4x16Decompress");
}

std::vector<std::byte> CramRans4x16Compress(std::span<const std::byte> data)
{
    if (std::empty(data)) {
        return {};
    }
    if (std::size(data) > std::numeric_limits<unsigned int>::max()) {
        throw std::runtime_error("CramRans4x16Compress: input size too large");
    }

    unsigned int outSize = 0;
    auto* out = rans_compress_4x16(NonConstBytes(data), static_cast<unsigned int>(std::size(data)),
                                   &outSize, 1);
    return CopyAndFree(out, outSize, "CramRans4x16Compress");
}

std::vector<std::byte> CramAdaptiveArithDecompress(std::span<const std::byte> data,
                                                   std::size_t rawSize)
{
    if (rawSize == 0) {
        return {};
    }
    if (std::size(data) > std::numeric_limits<unsigned int>::max() ||
        rawSize > std::numeric_limits<unsigned int>::max()) {
        throw std::runtime_error("CramAdaptiveArithDecompress: input/output size too large");
    }

    unsigned int outSize = 0;
    auto* out =
        arith_uncompress(NonConstBytes(data), static_cast<unsigned int>(std::size(data)), &outSize);
    if (outSize != rawSize) {
        std::free(out);
        throw std::runtime_error("CramAdaptiveArithDecompress: decompressed size mismatch");
    }
    return CopyAndFree(out, outSize, "CramAdaptiveArithDecompress");
}

std::vector<std::byte> CramAdaptiveArithCompress(std::span<const std::byte> data)
{
    if (std::empty(data)) {
        return {};
    }
    if (std::size(data) > std::numeric_limits<unsigned int>::max()) {
        throw std::runtime_error("CramAdaptiveArithCompress: input size too large");
    }

    unsigned int outSize = 0;
    auto* out = arith_compress(NonConstBytes(data), static_cast<unsigned int>(std::size(data)),
                               &outSize, 1);
    return CopyAndFree(out, outSize, "CramAdaptiveArithCompress");
}

std::vector<std::byte> CramFqzcompDecompress(std::span<const std::byte> data, std::size_t rawSize)
{
    if (rawSize == 0) {
        return {};
    }

    std::size_t outSize = 0;
    auto* out = fqz_decompress(reinterpret_cast<char*>(NonConstBytes(data)), std::size(data),
                               &outSize, nullptr, 0);
    if (outSize != rawSize) {
        std::free(out);
        throw std::runtime_error("CramFqzcompDecompress: decompressed size mismatch");
    }
    return CopyAndFree(out, outSize, "CramFqzcompDecompress");
}

std::vector<std::byte> CramFqzcompCompressSingleRecord(std::span<const std::byte> data)
{
    if (std::empty(data)) {
        return {};
    }
    if (std::size(data) > std::numeric_limits<std::uint32_t>::max()) {
        throw std::runtime_error("CramFqzcompCompressSingleRecord: input size too large");
    }

    std::uint32_t readLen = static_cast<std::uint32_t>(std::size(data));
    std::uint32_t flags = 0;
    return CramFqzcompCompress(data, std::span<const std::uint32_t>{&readLen, 1},
                               std::span<const std::uint32_t>{&flags, 1});
}

std::vector<std::byte> CramNameTokeniserDecompress(std::span<const std::byte> data,
                                                   std::size_t rawSize)
{
    if (rawSize == 0) {
        return {};
    }
    if (std::size(data) > std::numeric_limits<std::uint32_t>::max()) {
        throw std::runtime_error("CramNameTokeniserDecompress: input size too large");
    }

    std::uint32_t outSize = 0;
    auto* out = tok3_decode_names(NonConstBytes(data), static_cast<std::uint32_t>(std::size(data)),
                                  &outSize);
    if (outSize != rawSize) {
        std::free(out);
        throw std::runtime_error("CramNameTokeniserDecompress: decompressed size mismatch");
    }
    return CopyAndFree(out, outSize, "CramNameTokeniserDecompress");
}

std::vector<std::byte> CramNameTokeniserCompress(std::span<const std::byte> data)
{
    if (std::empty(data)) {
        return {};
    }
    if (std::size(data) > std::numeric_limits<int>::max()) {
        throw std::runtime_error("CramNameTokeniserCompress: input size too large");
    }

    int outLen = 0;
    int lastStart = 0;
    auto* out = tok3_encode_names(reinterpret_cast<char*>(NonConstBytes(data)),
                                  static_cast<int>(std::size(data)), 1, 0, &outLen, &lastStart);
    if (lastStart != static_cast<int>(std::size(data))) {
        std::free(out);
        throw std::runtime_error(
            "CramNameTokeniserCompress: input must contain complete separator-terminated names");
    }
    return CopyAndFree(out, static_cast<std::size_t>(outLen), "CramNameTokeniserCompress");
}

}  // namespace

std::vector<std::byte> CramFqzcompCompress(std::span<const std::byte> data,
                                           std::span<const std::uint32_t> recordLengths,
                                           std::span<const std::uint32_t> recordFlags)
{
    if (std::size(recordLengths) != std::size(recordFlags)) {
        throw std::runtime_error(
            "CramFqzcompCompress: recordLengths and recordFlags must have equal size");
    }
    if (std::size(recordLengths) > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::runtime_error("CramFqzcompCompress: record metadata size too large");
    }
    if (std::size(data) > std::numeric_limits<std::uint32_t>::max()) {
        throw std::runtime_error("CramFqzcompCompress: input size too large");
    }

    std::uint64_t totalLength = 0;
    for (const std::uint32_t length : recordLengths) {
        totalLength += length;
    }
    if (totalLength != std::size(data)) {
        throw std::runtime_error("CramFqzcompCompress: record lengths do not match input size");
    }
    if (std::empty(data)) {
        return {};
    }

    fqz_slice slice{.num_records = static_cast<int>(std::size(recordLengths)),
                    .len = const_cast<std::uint32_t*>(recordLengths.data()),
                    .flags = const_cast<std::uint32_t*>(recordFlags.data())};

    std::size_t outSize = 0;
    auto* out = fqz_compress(3 << 8, &slice, reinterpret_cast<char*>(NonConstBytes(data)),
                             std::size(data), &outSize, 0, nullptr);
    return CopyAndFree(out, outSize, "CramFqzcompCompress");
}

// ---------------------------------------------------------------------------
// Built-in gzip compression via libdeflate
// ---------------------------------------------------------------------------

std::vector<std::byte> CramGzipDecompress(std::span<const std::byte> data, std::size_t rawSize)
{
    return CramGzipDecompress(data, rawSize, nullptr);
}

std::vector<std::byte> CramGzipDecompress(std::span<const std::byte> data, std::size_t rawSize,
                                          libdeflate_decompressor* decompressor)
{
    LibdeflateDecompressorPtr ownedDecompressor{};
    auto* activeDecompressor = decompressor;
    if (activeDecompressor == nullptr) {
        ownedDecompressor.reset(libdeflate_alloc_decompressor());
        activeDecompressor = ownedDecompressor.get();
    }
    if (activeDecompressor == nullptr) {
        throw std::runtime_error("CramGzipDecompress: failed to allocate decompressor");
    }

    static std::once_flag libdeflateDispatchInitialized;
    std::call_once(libdeflateDispatchInitialized, [activeDecompressor, data, rawSize]() {
        // libdeflate lazily initializes a global decompression dispatch function.
        // Force that initialization once on a single thread before parallel decoding.
        std::vector<std::byte> warmupOutput(rawSize);
        std::size_t warmupOut = 0;
        (void)libdeflate_gzip_decompress(activeDecompressor, data.data(), std::size(data),
                                         warmupOutput.data(), rawSize, &warmupOut);
    });

    std::vector<std::byte> output(rawSize);
    std::size_t actualOut = 0;
    const auto result = libdeflate_gzip_decompress(activeDecompressor, data.data(), std::size(data),
                                                   output.data(), rawSize, &actualOut);

    if (result == LIBDEFLATE_SUCCESS && actualOut == rawSize) {
        return output;
    }

    // Try raw deflate if gzip fails (some CRAM blocks use raw deflate)
    const auto result2 = libdeflate_deflate_decompress(
        activeDecompressor, data.data(), std::size(data), output.data(), rawSize, &actualOut);

    if (result2 == LIBDEFLATE_SUCCESS && actualOut == rawSize) {
        return output;
    }

    throw std::runtime_error("CramGzipDecompress: decompression failed");
}

std::vector<std::byte> CramGzipCompress(std::span<const std::byte> data)
{
    return CramGzipCompress(data, nullptr);
}

std::vector<std::byte> CramGzipCompress(std::span<const std::byte> data,
                                        libdeflate_compressor* compressor)
{
    LibdeflateCompressorPtr ownedCompressor{};
    auto* activeCompressor = compressor;
    if (activeCompressor == nullptr) {
        ownedCompressor.reset(libdeflate_alloc_compressor(6));
        activeCompressor = ownedCompressor.get();
    }
    if (activeCompressor == nullptr) {
        throw std::runtime_error("CramGzipCompress: failed to allocate compressor");
    }

    const auto bound = libdeflate_gzip_compress_bound(activeCompressor, std::size(data));
    std::vector<std::byte> output(bound);

    const auto compressedSize = libdeflate_gzip_compress(
        activeCompressor, data.data(), std::size(data), output.data(), std::size(output));

    if (compressedSize == 0) {
        throw std::runtime_error("CramGzipCompress: compression failed");
    }

    output.resize(compressedSize);
    return output;
}

#ifdef PBSAMOA_HAVE_BZIP2
std::vector<std::byte> CramBzip2Decompress(std::span<const std::byte> data, std::size_t rawSize)
{
    if (rawSize > std::numeric_limits<unsigned int>::max() ||
        std::size(data) > std::numeric_limits<unsigned int>::max()) {
        throw std::runtime_error("CramBzip2Decompress: input/output size too large");
    }

    std::vector<std::byte> output(rawSize);
    unsigned int destLen = static_cast<unsigned int>(rawSize);
    const unsigned int srcLen = static_cast<unsigned int>(std::size(data));

    const auto rc = BZ2_bzBuffToBuffDecompress(
        reinterpret_cast<char*>(output.data()), &destLen,
        const_cast<char*>(reinterpret_cast<const char*>(data.data())), srcLen, 0, 0);
    if (rc != BZ_OK || destLen != rawSize) {
        throw std::runtime_error("CramBzip2Decompress: decompression failed");
    }
    return output;
}

std::vector<std::byte> CramBzip2Compress(std::span<const std::byte> data)
{
    if (std::empty(data)) {
        return {};
    }
    if (std::size(data) > std::numeric_limits<unsigned int>::max()) {
        throw std::runtime_error("CramBzip2Compress: input size too large");
    }
    const auto srcLen = static_cast<unsigned int>(std::size(data));
    // bzip2 manual: destination buffer should be source + 1% + 600.
    unsigned int destLen = srcLen + (srcLen / 100) + 601;
    std::vector<std::byte> output(destLen);

    const auto rc = BZ2_bzBuffToBuffCompress(
        reinterpret_cast<char*>(output.data()), &destLen,
        const_cast<char*>(reinterpret_cast<const char*>(data.data())), srcLen, 9, 0, 30);
    if (rc != BZ_OK) {
        throw std::runtime_error("CramBzip2Compress: compression failed");
    }
    output.resize(destLen);
    return output;
}
#endif

#ifdef PBSAMOA_HAVE_LZMA
std::vector<std::byte> CramLzmaDecompress(std::span<const std::byte> data, std::size_t rawSize)
{
    std::vector<std::byte> output(rawSize);
    std::size_t inPos = 0;
    std::size_t outPos = 0;
    std::uint64_t memlimit = UINT64_MAX;
    const std::uint32_t flags = 0;

    const auto rc = lzma_stream_buffer_decode(
        &memlimit, flags, nullptr, reinterpret_cast<const std::uint8_t*>(data.data()), &inPos,
        std::size(data), reinterpret_cast<std::uint8_t*>(output.data()), &outPos, rawSize);
    if (rc != LZMA_OK || inPos != std::size(data) || outPos != rawSize) {
        throw std::runtime_error("CramLzmaDecompress: decompression failed");
    }
    return output;
}

std::vector<std::byte> CramLzmaCompress(std::span<const std::byte> data)
{
    if (std::empty(data)) {
        return {};
    }
    const auto bound = lzma_stream_buffer_bound(std::size(data));
    std::vector<std::byte> output(bound);
    std::size_t outPos = 0;

    const auto rc = lzma_easy_buffer_encode(
        6, LZMA_CHECK_CRC32, nullptr, reinterpret_cast<const std::uint8_t*>(data.data()),
        std::size(data), reinterpret_cast<std::uint8_t*>(output.data()), &outPos,
        std::size(output));
    if (rc != LZMA_OK) {
        throw std::runtime_error("CramLzmaCompress: compression failed");
    }
    output.resize(outPos);
    return output;
}
#endif

// ---------------------------------------------------------------------------
// CompressionRegistry
// ---------------------------------------------------------------------------

CompressionRegistry::CompressionRegistry()
{
    // Register raw (no-op)
    methods_[0] = {
        .Compress = [](std::span<const std::byte> data) -> std::vector<std::byte> {
            return {data.begin(), data.end()};
        },
        .Decompress = [](std::span<const std::byte> data, std::size_t /*rawSize*/)
            -> std::vector<std::byte> { return {data.begin(), data.end()}; },
    };

    // Register gzip
    methods_[1] = {
        .Compress = [](std::span<const std::byte> data) { return CramGzipCompress(data); },
        .Decompress = [](std::span<const std::byte> data,
                         std::size_t rawSize) { return CramGzipDecompress(data, rawSize); },
    };

    methods_[4] = {
        .Compress = CramRans4x8Compress,
        .Decompress = CramRans4x8Decompress,
    };

    methods_[5] = {
        .Compress = CramRans4x16Compress,
        .Decompress = CramRans4x16Decompress,
    };

    methods_[6] = {
        .Compress = CramAdaptiveArithCompress,
        .Decompress = CramAdaptiveArithDecompress,
    };

    methods_[7] = {
        .Compress = CramFqzcompCompressSingleRecord,
        .Decompress = CramFqzcompDecompress,
    };

    methods_[8] = {
        .Compress = CramNameTokeniserCompress,
        .Decompress = CramNameTokeniserDecompress,
    };

#ifdef PBSAMOA_HAVE_BZIP2
    methods_[2] = {
        .Compress = CramBzip2Compress,
        .Decompress = CramBzip2Decompress,
    };
#endif

#ifdef PBSAMOA_HAVE_LZMA
    methods_[3] = {
        .Compress = CramLzmaCompress,
        .Decompress = CramLzmaDecompress,
    };
#endif
}

CompressionRegistry& CompressionRegistry::Instance()
{
    static CompressionRegistry instance;
    return instance;
}

void CompressionRegistry::Register(std::uint8_t methodId, CramCompressFn compress,
                                   CramDecompressFn decompress)
{
    if (methodId <= 8) {
        throw std::runtime_error("CompressionRegistry: cannot override built-in methods 0-8");
    }
    methods_[methodId] = {.Compress = std::move(compress), .Decompress = std::move(decompress)};
}

bool CompressionRegistry::HasMethod(std::uint8_t methodId) const
{
    return methods_.contains(methodId);
}

std::vector<std::byte> CompressionRegistry::Decompress(CramBlockMethod method,
                                                       std::span<const std::byte> data,
                                                       std::size_t rawSize) const
{
    const auto id = std::to_underlying(method);
    const auto it = methods_.find(id);
    if (it == methods_.end()) {
        throw std::runtime_error{std::format("CompressionRegistry: unsupported method {}", id)};
    }
    return it->second.Decompress(data, rawSize);
}

std::vector<std::byte> CompressionRegistry::Compress(CramBlockMethod method,
                                                     std::span<const std::byte> data) const
{
    const auto id = std::to_underlying(method);
    const auto it = methods_.find(id);
    if (it == methods_.end()) {
        throw std::runtime_error{std::format("CompressionRegistry: unsupported method {}", id)};
    }
    return it->second.Compress(data);
}

// ---------------------------------------------------------------------------
// Block-level helpers
// ---------------------------------------------------------------------------

void DecompressCramBlock(CramBlock& block) { DecompressCramBlock(block, nullptr); }

void DecompressCramBlock(CramBlock& block, libdeflate_decompressor* decompressor)
{
    if (block.RawSize == 0 || block.Method == CramBlockMethod::RAW) {
        return;  // already raw
    }
    if (block.Method == CramBlockMethod::GZIP) {
        block.Data =
            CramGzipDecompress(block.Data, static_cast<std::size_t>(block.RawSize), decompressor);
    } else {
        block.Data = CompressionRegistry::Instance().Decompress(
            block.Method, block.Data, static_cast<std::size_t>(block.RawSize));
    }
    block.CompressedSize = block.RawSize;
    block.Method = CramBlockMethod::RAW;
}

void CompressCramBlock(CramBlock& block, CramBlockMethod method)
{
    CompressCramBlock(block, method, nullptr);
}

void CompressCramBlock(CramBlock& block, CramBlockMethod method, libdeflate_compressor* compressor)
{
    if (method == CramBlockMethod::RAW) {
        block.Method = CramBlockMethod::RAW;
        block.CompressedSize = block.RawSize;
        return;
    }
    block.RawSize = static_cast<std::int32_t>(std::size(block.Data));
    if (method == CramBlockMethod::GZIP) {
        block.Data = CramGzipCompress(block.Data, compressor);
    } else {
        block.Data = CompressionRegistry::Instance().Compress(method, block.Data);
    }
    block.CompressedSize = static_cast<std::int32_t>(std::size(block.Data));
    block.Method = method;
}

}  // namespace Samoa
}  // namespace PacBio
