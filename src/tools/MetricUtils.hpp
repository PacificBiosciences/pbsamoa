#ifndef PBSAMOA_TOOLS_METRICUTILS_HPP
#define PBSAMOA_TOOLS_METRICUTILS_HPP

namespace PacBio {
namespace Samoa {
namespace Tools {

inline constexpr double BytesPerMiB{1024.0 * 1024.0};
inline constexpr double NsPerMs{1e6};

template <typename T>
double ToMiB(T bytes)
{
    return static_cast<double>(bytes) / BytesPerMiB;
}

template <typename T>
double ToMs(T nanoseconds)
{
    return static_cast<double>(nanoseconds) / NsPerMs;
}

template <typename T>
double ToRate(T count, double secs)
{
    return static_cast<double>(count) / secs;
}

template <typename T>
double RateOrZero(T count, double secs)
{
    if (secs <= 0) {
        return 0.0;
    }
    return ToRate(count, secs);
}

}  // namespace Tools
}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_TOOLS_METRICUTILS_HPP
