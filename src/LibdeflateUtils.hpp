#ifndef PBSAMOA_LIBDEFLATEUTILS_HPP
#define PBSAMOA_LIBDEFLATEUTILS_HPP

#include <libdeflate.h>

#include <memory>

namespace PacBio {
namespace Samoa {

struct LibdeflateDecompressorDeleter
{
    void operator()(libdeflate_decompressor* d) const noexcept { libdeflate_free_decompressor(d); }
};

struct LibdeflateCompressorDeleter
{
    void operator()(libdeflate_compressor* c) const noexcept { libdeflate_free_compressor(c); }
};

using LibdeflateDecompressorPtr =
    std::unique_ptr<libdeflate_decompressor, LibdeflateDecompressorDeleter>;
using LibdeflateCompressorPtr = std::unique_ptr<libdeflate_compressor, LibdeflateCompressorDeleter>;

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_LIBDEFLATEUTILS_HPP
