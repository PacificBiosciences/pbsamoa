#ifndef PBSAMOA_CORE_GENOMICINTERVAL_HPP
#define PBSAMOA_CORE_GENOMICINTERVAL_HPP

#include <string>
#include <string_view>

#include <cstdint>

namespace PacBio {
namespace Samoa {

/// \brief Half-open genomic interval [start, stop) on a named reference.
class GenomicInterval
{
public:
    GenomicInterval();

    GenomicInterval(std::string name, std::int32_t start, std::int32_t stop);

    std::string_view Name() const;

    std::int32_t Start() const;

    std::int32_t Stop() const;

    bool Empty() const;

    GenomicInterval& Name(std::string name);

    GenomicInterval& Start(std::int32_t start);

    GenomicInterval& Stop(std::int32_t stop);

private:
    static void Validate(std::int32_t start, std::int32_t stop);

    std::string name_;
    std::int32_t start_{0};
    std::int32_t stop_{0};
};

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_CORE_GENOMICINTERVAL_HPP
