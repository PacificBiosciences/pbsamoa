#include "ZmwUtils.hpp"

#include <charconv>
#include <system_error>

namespace PacBio {
namespace Samoa {
namespace {

std::string_view PrefixBeforeSlash(std::string_view text) { return text.substr(0, text.find('/')); }

template <typename Integer>
Integer ParseOrZero(std::string_view text, int base = 10)
{
    Integer value{0};
    const std::from_chars_result parseResult{
        std::from_chars(std::data(text), std::data(text) + std::size(text), value, base)};
    if (parseResult.ec != std::errc{}) {
        return 0;
    }
    return value;
}

}  // namespace

std::int32_t ParseZmwFromName(std::string_view name)
{
    const std::size_t firstSlash{name.find('/')};
    if (firstSlash == std::string_view::npos) {
        return 0;
    }
    return ParseOrZero<std::int32_t>(PrefixBeforeSlash(name.substr(firstSlash + 1)));
}

std::string_view MovieZmwPrefix(std::string_view name)
{
    const std::size_t firstSlash{name.find('/')};
    if (firstSlash == std::string_view::npos) {
        return name;
    }
    const std::string_view zmwField{PrefixBeforeSlash(name.substr(firstSlash + 1))};
    if ((firstSlash + 1 + std::size(zmwField)) == std::size(name)) {
        return name;
    }
    return name.substr(0, firstSlash + 1 + std::size(zmwField));
}

std::string_view ReadGroupBaseId(std::string_view readGroupId)
{
    return PrefixBeforeSlash(readGroupId);
}

std::int32_t ParseReadGroupId(std::string_view readGroupId)
{
    const std::string_view baseId{ReadGroupBaseId(readGroupId)};
    if (baseId.empty()) {
        return 0;
    }
    return static_cast<std::int32_t>(ParseOrZero<std::uint32_t>(baseId, 16));
}

}  // namespace Samoa
}  // namespace PacBio
