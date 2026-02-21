#ifndef PBSAMOA_ZMWUTILS_HPP
#define PBSAMOA_ZMWUTILS_HPP

#include <charconv>
#include <string_view>
#include <system_error>

#include <cstdint>

namespace PacBio {
namespace Samoa {

/// \brief Parse ZMW hole number from PacBio read name format "movie/<zmw>/...".
/// \returns hole number, or 0 on parse failure.
inline std::int32_t ParseZmwFromName(std::string_view name)
{
    const std::size_t firstSlash{name.find('/')};
    if (firstSlash == std::string_view::npos) {
        return 0;
    }
    const std::string_view rest{name.substr(firstSlash + 1)};
    const std::size_t secondSlash{rest.find('/')};
    const std::string_view zmwStr{rest.substr(0, secondSlash)};

    std::int32_t zmw{0};
    const auto [ptr, ec] =
        std::from_chars(std::data(zmwStr), std::data(zmwStr) + std::size(zmwStr), zmw);
    if (ec != std::errc{}) {
        return 0;
    }
    return zmw;
}

/// \brief Extract the "movie/zmw" prefix from a PacBio read name.
/// \returns "movie/zmw" substring, or the full name if no slashes found.
inline std::string_view MovieZmwPrefix(std::string_view name)
{
    const std::size_t firstSlash{name.find('/')};
    if (firstSlash == std::string_view::npos) {
        return name;
    }
    const std::string_view rest{name.substr(firstSlash + 1)};
    const std::size_t secondSlash{rest.find('/')};
    if (secondSlash == std::string_view::npos) {
        return name;
    }
    return name.substr(0, firstSlash + 1 + secondSlash);
}

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_ZMWUTILS_HPP
