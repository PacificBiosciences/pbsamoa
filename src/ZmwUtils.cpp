#include "ZmwUtils.hpp"

#include <pbsamoa/core/SamHeader.hpp>
#include <pbsamoa/core/Tags.hpp>

#include <charconv>
#include <cstddef>
#include <system_error>
#include <variant>

namespace PacBio {
namespace Samoa {
namespace {

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

std::int32_t ReadGroupIdFromTags(const TagMap& tags)
{
    const TagValue* rgValue{tags.Get(RG_TAG)};
    if (const auto* rgText{std::get_if<std::string>(rgValue)}) {
        return ParseReadGroupId(*rgText);
    }
    return 0;
}

}  // namespace

std::int32_t ParseZmwFromName(std::string_view name)
{
    const std::size_t firstSlash{name.find('/')};
    if (firstSlash == std::string_view::npos) {
        return 0;
    }

    const std::size_t secondSlash{name.find('/', firstSlash + 1)};
    const std::size_t zmwStart{firstSlash + 1};
    std::string_view zmwField{};
    if (secondSlash == std::string_view::npos) {
        zmwField = name.substr(zmwStart);
    } else {
        zmwField = name.substr(zmwStart, secondSlash - zmwStart);
    }
    if (zmwField.empty()) {
        return 0;
    }
    return ParseOrZero<std::int32_t>(zmwField);
}

std::int32_t ParseReadGroupId(std::string_view readGroupId)
{
    const std::string_view baseId{ReadGroupBaseId(readGroupId)};
    if (baseId.empty()) {
        return 0;
    }
    return static_cast<std::int32_t>(ParseOrZero<std::uint32_t>(baseId, 16));
}

ZmwIdentity ParseZmwIdentity(std::string_view name, const TagMap& tags)
{
    return ZmwIdentity{ReadGroupIdFromTags(tags), ParseZmwFromName(name)};
}

}  // namespace Samoa
}  // namespace PacBio
