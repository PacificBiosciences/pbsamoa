#include <pbsamoa/core/GenomicInterval.hpp>

#include <stdexcept>
#include <utility>

namespace PacBio {
namespace Samoa {

GenomicInterval::GenomicInterval() = default;

GenomicInterval::GenomicInterval(std::string name, std::int32_t start, std::int32_t stop)
    : name_{std::move(name)}, start_{start}, stop_{stop}
{
    Validate(start_, stop_);
}

std::string_view GenomicInterval::Name() const { return name_; }

std::int32_t GenomicInterval::Start() const { return start_; }

std::int32_t GenomicInterval::Stop() const { return stop_; }

bool GenomicInterval::Empty() const { return start_ == stop_; }

GenomicInterval& GenomicInterval::Name(std::string name)
{
    name_ = std::move(name);
    return *this;
}

GenomicInterval& GenomicInterval::Start(std::int32_t start)
{
    Validate(start, stop_);
    start_ = start;
    return *this;
}

GenomicInterval& GenomicInterval::Stop(std::int32_t stop)
{
    Validate(start_, stop);
    stop_ = stop;
    return *this;
}

void GenomicInterval::Validate(std::int32_t start, std::int32_t stop)
{
    if (start < 0) {
        throw std::invalid_argument{"GenomicInterval: start must be >= 0"};
    }
    if (stop < start) {
        throw std::invalid_argument{"GenomicInterval: stop must be >= start"};
    }
}

}  // namespace Samoa
}  // namespace PacBio
