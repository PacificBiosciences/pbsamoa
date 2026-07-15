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
#include <array>
#include <format>
#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string_view>
#include <utility>

#include <cstdlib>

namespace PacBio {
namespace Samoa {

std::vector<std::byte> CramFqzcompCompress(std::span<const std::byte> data,
                                           std::span<const std::uint32_t> recordLengths,
                                           std::span<const std::uint32_t> recordFlags);

namespace {

std::vector<std::byte> CopyRawBlock(std::span<const std::byte> data)
{
    return {std::begin(data), std::end(data)};
}

std::vector<std::byte> DecompressRawBlock(std::span<const std::byte> data,
                                          std::size_t /* rawSize */)
{
    return CopyRawBlock(data);
}

template <typename T>
T CheckedSizeCast(const std::size_t size, const std::string_view context,
                  const std::string_view reason)
{
    if (size > static_cast<std::size_t>(std::numeric_limits<T>::max())) {
        throw std::runtime_error(std::format("{}: {}", context, reason));
    }
    return static_cast<T>(size);
}

std::vector<std::byte> CopyAndFree(unsigned char* ptr, std::size_t size, std::string_view context)
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

unsigned char* NonConstBytes(std::span<const std::byte> data)
{
    return reinterpret_cast<unsigned char*>(const_cast<std::byte*>(data.data()));
}

template <typename SizeType, typename DecompressFn>
std::vector<std::byte> CheckedDecompress(std::span<const std::byte> data, std::size_t rawSize,
                                         std::string_view context, DecompressFn decompressFn)
{
    if (rawSize == 0) {
        return {};
    }

    const SizeType inputSize{
        CheckedSizeCast<SizeType>(std::size(data), context, "input/output size too large")};
    const SizeType expectedOutputSize{
        CheckedSizeCast<SizeType>(rawSize, context, "input/output size too large")};

    SizeType outSize{0};
    auto* out = decompressFn(NonConstBytes(data), inputSize, &outSize);
    if (outSize != expectedOutputSize) {
        std::free(out);
        throw std::runtime_error(std::format("{}: decompressed size mismatch", context));
    }
    return CopyAndFree(out, static_cast<std::size_t>(outSize), context);
}

template <typename SizeType, typename CompressFn>
std::vector<std::byte> CheckedCompress(std::span<const std::byte> data, std::string_view context,
                                       CompressFn compressFn)
{
    if (std::empty(data)) {
        return {};
    }

    const SizeType inputSize{
        CheckedSizeCast<SizeType>(std::size(data), context, "input size too large")};

    SizeType outSize{0};
    auto* out = compressFn(NonConstBytes(data), inputSize, &outSize);
    return CopyAndFree(out, static_cast<std::size_t>(outSize), context);
}

unsigned char* Rans4x8CompressBytes(unsigned char* input, unsigned int inputSize,
                                    unsigned int* outSize)
{
    return rans_compress(input, inputSize, outSize, 1);
}

unsigned char* Rans4x16CompressBytes(unsigned char* input, unsigned int inputSize,
                                     unsigned int* outSize)
{
    return rans_compress_4x16(input, inputSize, outSize, 1);
}

unsigned char* AdaptiveArithCompressBytes(unsigned char* input, unsigned int inputSize,
                                          unsigned int* outSize)
{
    return arith_compress(input, inputSize, outSize, 1);
}

std::vector<std::byte> CramRans4x8Decompress(std::span<const std::byte> data, std::size_t rawSize)
{
    return CheckedDecompress<unsigned int>(data, rawSize, "CramRans4x8Decompress", rans_uncompress);
}

std::vector<std::byte> CramRans4x8Compress(std::span<const std::byte> data)
{
    return CheckedCompress<unsigned int>(data, "CramRans4x8Compress", Rans4x8CompressBytes);
}

std::vector<std::byte> CramRans4x16Decompress(std::span<const std::byte> data, std::size_t rawSize)
{
    return CheckedDecompress<unsigned int>(data, rawSize, "CramRans4x16Decompress",
                                           rans_uncompress_4x16);
}

std::vector<std::byte> CramRans4x16Compress(std::span<const std::byte> data)
{
    return CheckedCompress<unsigned int>(data, "CramRans4x16Compress", Rans4x16CompressBytes);
}

std::vector<std::byte> CramAdaptiveArithDecompress(std::span<const std::byte> data,
                                                   std::size_t rawSize)
{
    return CheckedDecompress<unsigned int>(data, rawSize, "CramAdaptiveArithDecompress",
                                           arith_uncompress);
}

std::vector<std::byte> CramAdaptiveArithCompress(std::span<const std::byte> data)
{
    return CheckedCompress<unsigned int>(data, "CramAdaptiveArithCompress",
                                         AdaptiveArithCompressBytes);
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
    return CopyAndFree(reinterpret_cast<unsigned char*>(out), outSize, "CramFqzcompDecompress");
}

std::vector<std::byte> CramFqzcompCompressSingleRecord(std::span<const std::byte> data)
{
    if (std::empty(data)) {
        return {};
    }
    const std::uint32_t readLen{CheckedSizeCast<std::uint32_t>(
        std::size(data), "CramFqzcompCompressSingleRecord", "input size too large")};
    const std::array<std::uint32_t, 1> recordLengths{readLen};
    const std::array<std::uint32_t, 1> recordFlags{0};
    return CramFqzcompCompress(data, recordLengths, recordFlags);
}

std::vector<std::byte> CramNameTokeniserDecompress(std::span<const std::byte> data,
                                                   std::size_t rawSize)
{
    if (rawSize == 0) {
        return {};
    }
    const std::uint32_t inputSize{CheckedSizeCast<std::uint32_t>(
        std::size(data), "CramNameTokeniserDecompress", "input size too large")};

    std::uint32_t outSize{0};
    auto* out = tok3_decode_names(NonConstBytes(data), inputSize, &outSize);
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
    const int inputSize{
        CheckedSizeCast<int>(std::size(data), "CramNameTokeniserCompress", "input size too large")};

    int outLen = 0;
    int lastStart = 0;
    auto* out = tok3_encode_names(reinterpret_cast<char*>(NonConstBytes(data)), inputSize, 1, 0,
                                  &outLen, &lastStart);
    if (lastStart != inputSize) {
        std::free(out);
        throw std::runtime_error(
            "CramNameTokeniserCompress: input must contain "
            "complete separator-terminated names");
    }
    return CopyAndFree(out, static_cast<std::size_t>(outLen), "CramNameTokeniserCompress");
}

void WarmLibdeflateDispatch(libdeflate_decompressor* decompressor, std::span<const std::byte> data,
                            std::size_t rawSize)
{
    // libdeflate lazily initializes a global decompression dispatch
    // function. Force that initialization once on a single thread before
    // parallel decoding.
    std::vector<std::byte> warmupOutput(rawSize);
    std::size_t warmupOut = 0;
    (void)libdeflate_gzip_decompress(decompressor, data.data(), std::size(data),
                                     warmupOutput.data(), rawSize, &warmupOut);
}

}  // namespace

std::vector<std::byte> CramFqzcompCompress(std::span<const std::byte> data,
                                           std::span<const std::uint32_t> recordLengths,
                                           std::span<const std::uint32_t> recordFlags)
{
    if (std::size(recordLengths) != std::size(recordFlags)) {
        throw std::runtime_error(
            "CramFqzcompCompress: recordLengths and "
            "recordFlags must have equal size");
    }
    const int numRecords{CheckedSizeCast<int>(std::size(recordLengths), "CramFqzcompCompress",
                                              "record metadata size too large")};
    const std::uint32_t inputSize{CheckedSizeCast<std::uint32_t>(
        std::size(data), "CramFqzcompCompress", "input size too large")};

    const std::uint64_t totalLength =
        std::ranges::fold_left(recordLengths, std::uint64_t{0}, std::plus{});
    if (totalLength != std::size(data)) {
        throw std::runtime_error("CramFqzcompCompress: record lengths do not match input size");
    }
    if (std::empty(data)) {
        return {};
    }

    fqz_slice slice{.num_records = numRecords,
                    .len = const_cast<std::uint32_t*>(recordLengths.data()),
                    .flags = const_cast<std::uint32_t*>(recordFlags.data())};

    std::size_t outSize = 0;
    auto* out = fqz_compress(3 << 8, &slice, reinterpret_cast<char*>(NonConstBytes(data)),
                             inputSize, &outSize, 0, nullptr);
    return CopyAndFree(reinterpret_cast<unsigned char*>(out), outSize, "CramFqzcompCompress");
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
    if (!activeDecompressor) {
        ownedDecompressor.reset(libdeflate_alloc_decompressor());
        activeDecompressor = ownedDecompressor.get();
    }
    if (!activeDecompressor) {
        throw std::runtime_error("CramGzipDecompress: failed to allocate decompressor");
    }

    static std::once_flag libdeflateDispatchInitialized;
    std::call_once(libdeflateDispatchInitialized, WarmLibdeflateDispatch, activeDecompressor, data,
                   rawSize);

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
    return CramGzipCompress(data, nullptr, std::nullopt);
}

std::vector<std::byte> CramGzipCompress(std::span<const std::byte> data, const int compressionLevel)
{
    return CramGzipCompress(data, nullptr, compressionLevel);
}

std::vector<std::byte> CramGzipCompress(std::span<const std::byte> data,
                                        libdeflate_compressor* compressor,
                                        std::optional<int> compressionLevel)
{
    LibdeflateCompressorPtr ownedCompressor{};
    auto* activeCompressor = compressor;
    if (!activeCompressor) {
        const int level = compressionLevel.value_or(DEFAULT_GZIP_COMPRESSION_LEVEL);
        if (level < 0 || level > 12) {
            throw std::runtime_error(
                std::format("CramGzipCompress: invalid gzip compression level {}", level));
        }
        ownedCompressor.reset(libdeflate_alloc_compressor(level));
        activeCompressor = ownedCompressor.get();
    }
    if (!activeCompressor) {
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
    std::vector<std::byte> output(rawSize);
    unsigned int destLen = CheckedSizeCast<unsigned int>(rawSize, "CramBzip2Decompress",
                                                         "input/output size too large");
    const unsigned int srcLen = CheckedSizeCast<unsigned int>(
        std::size(data), "CramBzip2Decompress", "input/output size too large");

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
    const auto srcLen =
        CheckedSizeCast<unsigned int>(std::size(data), "CramBzip2Compress", "input size too large");
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
        .Compress = CopyRawBlock,
        .Decompress = DecompressRawBlock,
    };

    // Register gzip
    methods_[1] = {
        .Compress =
            static_cast<std::vector<std::byte> (*)(std::span<const std::byte>)>(&CramGzipCompress),
        .Decompress =
            static_cast<std::vector<std::byte> (*)(std::span<const std::byte>, std::size_t)>(
                &CramGzipDecompress),
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

const CompressionRegistry::MethodEntry& CompressionRegistry::LookupMethod(
    const std::uint8_t methodId) const
{
    const auto it = methods_.find(methodId);
    if (it == methods_.end()) {
        throw std::runtime_error{
            std::format("CompressionRegistry: unsupported method {}", methodId)};
    }
    return it->second;
}

std::vector<std::byte> CompressionRegistry::Decompress(CramBlockMethod method,
                                                       std::span<const std::byte> data,
                                                       std::size_t rawSize) const
{
    return LookupMethod(std::to_underlying(method)).Decompress(data, rawSize);
}

std::vector<std::byte> CompressionRegistry::Compress(CramBlockMethod method,
                                                     std::span<const std::byte> data) const
{
    return LookupMethod(std::to_underlying(method)).Compress(data);
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
    switch (block.Method) {
        case CramBlockMethod::GZIP:
            block.Data = CramGzipDecompress(block.Data, static_cast<std::size_t>(block.RawSize),
                                            decompressor);
            break;
        default:
            block.Data = CompressionRegistry::Instance().Decompress(
                block.Method, block.Data, static_cast<std::size_t>(block.RawSize));
            break;
    }
    block.CompressedSize = block.RawSize;
    block.Method = CramBlockMethod::RAW;
}

void CompressCramBlock(CramBlock& block, CramBlockMethod method,
                       std::optional<int> compressionLevel)
{
    CompressCramBlock(block, method, nullptr, compressionLevel);
}

void CompressCramBlock(CramBlock& block, CramBlockMethod method, libdeflate_compressor* compressor,
                       std::optional<int> compressionLevel)
{
    if (method == CramBlockMethod::RAW) {
        block.Method = CramBlockMethod::RAW;
        block.CompressedSize = block.RawSize;
        return;
    }
    block.RawSize = static_cast<std::int32_t>(std::size(block.Data));
    switch (method) {
        case CramBlockMethod::GZIP:
            block.Data = CramGzipCompress(block.Data, compressor, compressionLevel);
            break;
        default:
            block.Data = CompressionRegistry::Instance().Compress(method, block.Data);
            break;
    }
    block.CompressedSize = static_cast<std::int32_t>(std::size(block.Data));
    block.Method = method;
}

}  // namespace Samoa
}  // namespace PacBio
