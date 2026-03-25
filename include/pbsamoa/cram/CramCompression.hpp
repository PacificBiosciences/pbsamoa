#ifndef PBSAMOA_CRAM_CRAMCOMPRESSION_HPP
#define PBSAMOA_CRAM_CRAMCOMPRESSION_HPP

#include <pbsamoa/cram/CramStructs.hpp>

#include <libdeflate.h>

#include <functional>
#include <optional>
#include <span>
#include <unordered_map>
#include <vector>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {

// ---------------------------------------------------------------------------
// Compression / Decompression function types
// ---------------------------------------------------------------------------

/// \brief Decompression function: (compressed_data, raw_size) -> decompressed
/// bytes. Must produce exactly rawSize bytes.
using CramDecompressFn =
    std::function<std::vector<std::byte>(std::span<const std::byte> data, std::size_t rawSize)>;

/// \brief Compression function: (raw_data) -> compressed bytes.
using CramCompressFn = std::function<std::vector<std::byte>(std::span<const std::byte> data)>;

// ---------------------------------------------------------------------------
// CompressionRegistry
// ---------------------------------------------------------------------------

/// \brief Registry for CRAM block compression methods.
///
/// CRAM v3.x block methods (0-8) are registered by default. Users can register
/// additional custom methods by ID.
///
/// Thread-safe for reads after initial registration. Not thread-safe for
/// concurrent registration and use.
class CompressionRegistry
{
public:
    /// \brief Get the global singleton registry.
    static CompressionRegistry& Instance();

    /// \brief Register a custom compression method.
    /// \param methodId CRAM block method ID (0-255). Built-in CRAM IDs (0-8)
    /// cannot be overridden.
    /// \param compress compression function
    /// \param decompress decompression function
    void Register(std::uint8_t methodId, CramCompressFn compress, CramDecompressFn decompress);

    /// \brief Check if a method is registered.
    bool HasMethod(std::uint8_t methodId) const;

    /// \brief Decompress a block's data using the registered method.
    std::vector<std::byte> Decompress(CramBlockMethod method, std::span<const std::byte> data,
                                      std::size_t rawSize) const;

    /// \brief Compress data using the specified method.
    std::vector<std::byte> Compress(CramBlockMethod method, std::span<const std::byte> data) const;

private:
    CompressionRegistry();

    struct MethodEntry
    {
        CramCompressFn Compress;
        CramDecompressFn Decompress;
    };

    const MethodEntry& LookupMethod(std::uint8_t methodId) const;

    std::unordered_map<std::uint8_t, MethodEntry> methods_;
};

// ---------------------------------------------------------------------------
// Built-in compression helpers
// ---------------------------------------------------------------------------

/// \brief Decompress gzip data to rawSize bytes using libdeflate.
std::vector<std::byte> CramGzipDecompress(std::span<const std::byte> data, std::size_t rawSize);

/// \brief Decompress gzip data with a reusable libdeflate decompressor context.
std::vector<std::byte> CramGzipDecompress(std::span<const std::byte> data, std::size_t rawSize,
                                          libdeflate_decompressor* decompressor);

/// \brief Compress data using gzip via libdeflate.
std::vector<std::byte> CramGzipCompress(std::span<const std::byte> data);

/// \brief Compress data using gzip via libdeflate with an explicit level.
std::vector<std::byte> CramGzipCompress(std::span<const std::byte> data, int compressionLevel);

/// \brief Compress gzip data with a reusable libdeflate compressor context.
std::vector<std::byte> CramGzipCompress(std::span<const std::byte> data,
                                        libdeflate_compressor* compressor,
                                        std::optional<int> compressionLevel = std::nullopt);

/// \brief Compress quality data using fqzcomp with explicit per-record
/// metadata.
std::vector<std::byte> CramFqzcompCompress(std::span<const std::byte> data,
                                           std::span<const std::uint32_t> recordLengths,
                                           std::span<const std::uint32_t> recordFlags);

/// \brief Decompress a CRAM block in-place (replaces Data with uncompressed
/// content).
void DecompressCramBlock(CramBlock& block);

/// \brief Decompress a CRAM block in-place with an optional reusable gzip
/// context.
void DecompressCramBlock(CramBlock& block, libdeflate_decompressor* decompressor);

/// \brief Compress a CRAM block in-place for writing (replaces Data with
/// compressed content).
void CompressCramBlock(CramBlock& block, CramBlockMethod method,
                       std::optional<int> compressionLevel = std::nullopt);

/// \brief Compress a CRAM block in-place with an optional reusable gzip
/// context.
void CompressCramBlock(CramBlock& block, CramBlockMethod method, libdeflate_compressor* compressor,
                       std::optional<int> compressionLevel = std::nullopt);

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_CRAM_CRAMCOMPRESSION_HPP
