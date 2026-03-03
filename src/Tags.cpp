#include <pbsamoa/core/Tags.hpp>

#include "BinaryUtils.hpp"
#include "CramInternal.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <charconv>
#include <format>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

#include <cstdio>

namespace PacBio {
namespace Samoa {

namespace {

void WriteU8(std::vector<std::byte>& out, std::uint8_t v)
{
    out.push_back(static_cast<std::byte>(v));
}

void WriteI8(std::vector<std::byte>& out, std::int8_t v)
{
    out.push_back(static_cast<std::byte>(v));
}

void AppendU16LE(std::vector<std::byte>& out, std::uint16_t v)
{
    out.push_back(static_cast<std::byte>(v & 0xFFU));
    out.push_back(static_cast<std::byte>((v >> 8U) & 0xFFU));
}

void AppendU32LE(std::vector<std::byte>& out, std::uint32_t v)
{
    out.push_back(static_cast<std::byte>(v & 0xFFU));
    out.push_back(static_cast<std::byte>((v >> 8U) & 0xFFU));
    out.push_back(static_cast<std::byte>((v >> 16U) & 0xFFU));
    out.push_back(static_cast<std::byte>((v >> 24U) & 0xFFU));
}

void AppendI16LE(std::vector<std::byte>& out, std::int16_t v)
{
    AppendU16LE(out, std::bit_cast<std::uint16_t>(v));
}

void AppendI32LE(std::vector<std::byte>& out, std::int32_t v)
{
    AppendU32LE(out, std::bit_cast<std::uint32_t>(v));
}

void AppendF32LE(std::vector<std::byte>& out, float v)
{
    AppendU32LE(out, std::bit_cast<std::uint32_t>(v));
}

std::int16_t ReadI16LE(const std::byte* data)
{
    return std::bit_cast<std::int16_t>(ReadU16LE(data));
}

float ReadF32LE(const std::byte* data) { return std::bit_cast<float>(ReadU32LE(data)); }

char SmallestIntType(std::int64_t v)
{
    if ((v >= 0) && (v <= 255)) {
        return 'C';
    }
    if ((v >= -128) && (v <= 127)) {
        return 'c';
    }
    if ((v >= 0) && (v <= 65535)) {
        return 'S';
    }
    if ((v >= -32768) && (v <= 32767)) {
        return 's';
    }
    if ((v >= 0) && (v <= 4294967295LL)) {
        return 'I';
    }
    return 'i';
}

bool IsHexDigit(char c)
{
    return ((c >= '0') && (c <= '9')) || ((c >= 'a') && (c <= 'f')) || ((c >= 'A') && (c <= 'F'));
}

void RequirePayloadSize(char type, std::span<const std::byte> payload, std::size_t expected)
{
    if (std::size(payload) != expected) {
        throw std::runtime_error{std::format("Tags: type '{}' requires {} payload byte(s), got {}",
                                             type, expected, std::size(payload))};
    }
}

struct BamSerializeVisitor
{
    std::vector<std::byte>& result;

    void operator()(char v) const
    {
        result.push_back(static_cast<std::byte>('A'));
        result.push_back(static_cast<std::byte>(v));
    }

    void operator()(std::int64_t v) const
    {
        const char bamType{SmallestIntType(v)};
        result.push_back(static_cast<std::byte>(bamType));
        switch (bamType) {
            case 'c':
                WriteI8(result, v);
                break;
            case 'C':
                WriteU8(result, v);
                break;
            case 's': {
                const std::int16_t sv = v;
                AppendI16LE(result, sv);
                break;
            }
            case 'S': {
                const std::uint16_t sv = v;
                AppendU16LE(result, sv);
                break;
            }
            case 'i': {
                const std::int32_t sv = v;
                AppendI32LE(result, sv);
                break;
            }
            case 'I': {
                const std::uint32_t sv = v;
                AppendU32LE(result, sv);
                break;
            }
            default:
                break;
        }
    }

    void operator()(float v) const
    {
        result.push_back(static_cast<std::byte>('f'));
        AppendF32LE(result, v);
    }

    void operator()(std::string_view v) const
    {
        result.push_back(static_cast<std::byte>('Z'));
        for (const char ch : v) {
            result.push_back(static_cast<std::byte>(ch));
        }
        result.push_back(std::byte{0});  // NUL terminator
    }

    void operator()(const HexString& v) const
    {
        result.push_back(static_cast<std::byte>('H'));
        for (const char ch : v.value) {
            result.push_back(static_cast<std::byte>(ch));
        }
        result.push_back(std::byte{0});  // NUL terminator
    }

    void operator()(const TagArray& v) const
    {
        result.push_back(static_cast<std::byte>('B'));
        result.push_back(static_cast<std::byte>(v.ElementType()));
        const std::uint32_t count{v.Count()};
        AppendU32LE(result, count);
        const std::span<const std::byte> data{v.Data()};
        result.insert(std::ranges::end(result), std::ranges::begin(data), std::ranges::end(data));
    }
};

struct SamSerializeVisitor
{
    std::string& result;

    void operator()(char v) const
    {
        result += "A:";
        result += v;
    }

    void operator()(std::int64_t v) const
    {
        result += "i:";
        std::format_to(std::back_inserter(result), "{}", v);
    }

    void operator()(float v) const
    {
        result += "f:";
        std::format_to(std::back_inserter(result), "{:g}", v);
    }

    void operator()(std::string_view v) const
    {
        result += "Z:";
        result += v;
    }

    void operator()(const HexString& v) const
    {
        result += "H:";
        result += v.value;
    }

    void operator()(const TagArray& v) const
    {
        result += "B:";
        result += v.ElementType();

        const std::size_t count{v.Count()};
        const std::span<const std::byte> data{v.Data()};
        const char elemType{v.ElementType()};

        for (std::size_t i{0}; i < count; ++i) {
            result += ',';
            if (elemType == 'f') {
                const float fv{ReadF32LE(std::data(data) + (i * 4))};
                std::format_to(std::back_inserter(result), "{:g}", fv);
            } else if (elemType == 'c') {
                const std::int8_t iv{static_cast<std::int8_t>(data[i])};
                std::format_to(std::back_inserter(result), "{}", iv);
            } else if (elemType == 'C') {
                const std::uint8_t iv{static_cast<std::uint8_t>(data[i])};
                std::format_to(std::back_inserter(result), "{}", iv);
            } else if (elemType == 's') {
                const std::int16_t iv{ReadI16LE(std::data(data) + (i * 2))};
                std::format_to(std::back_inserter(result), "{}", iv);
            } else if (elemType == 'S') {
                const std::uint16_t iv{ReadU16LE(std::data(data) + (i * 2))};
                std::format_to(std::back_inserter(result), "{}", iv);
            } else if (elemType == 'i') {
                const std::int32_t iv{ReadI32LE(std::data(data) + (i * 4))};
                std::format_to(std::back_inserter(result), "{}", iv);
            } else if (elemType == 'I') {
                const std::uint32_t iv{ReadU32LE(std::data(data) + (i * 4))};
                std::format_to(std::back_inserter(result), "{}", iv);
            }
        }
    }
};

}  // namespace

TagValue DecodeTagValueFromBamPayload(char type, std::span<const std::byte> payload)
{
    switch (type) {
        case 'A':
            RequirePayloadSize(type, payload, 1);
            return TagValue{static_cast<char>(std::to_integer<std::uint8_t>(payload.front()))};
        case 'c':
            RequirePayloadSize(type, payload, 1);
            return TagValue{std::int64_t{
                static_cast<std::int8_t>(std::to_integer<std::uint8_t>(payload.front()))}};
        case 'C':
            RequirePayloadSize(type, payload, 1);
            return TagValue{std::int64_t{std::to_integer<std::uint8_t>(payload.front())}};
        case 's':
            RequirePayloadSize(type, payload, 2);
            return TagValue{std::int64_t{ReadI16LE(payload.data())}};
        case 'S':
            RequirePayloadSize(type, payload, 2);
            return TagValue{std::int64_t{ReadU16LE(payload.data())}};
        case 'i':
            RequirePayloadSize(type, payload, 4);
            return TagValue{std::int64_t{ReadI32LE(payload.data())}};
        case 'I':
            RequirePayloadSize(type, payload, 4);
            return TagValue{std::int64_t{ReadU32LE(payload.data())}};
        case 'f':
            RequirePayloadSize(type, payload, 4);
            return TagValue{ReadF32LE(payload.data())};
        case 'Z': {
            std::size_t textLen = std::size(payload);
            if (textLen > 0 && payload[textLen - 1] == std::byte{0}) {
                --textLen;
            }
            return TagValue{std::string{reinterpret_cast<const char*>(payload.data()), textLen}};
        }
        case 'H': {
            std::size_t textLen = std::size(payload);
            if (textLen > 0 && payload[textLen - 1] == std::byte{0}) {
                --textLen;
            }
            return TagValue{
                HexString{std::string{reinterpret_cast<const char*>(payload.data()), textLen}}};
        }
        case 'B': {
            if (std::size(payload) < 5) {
                throw std::runtime_error{"Tags: type 'B' payload too short"};
            }
            const char elemType = static_cast<char>(std::to_integer<std::uint8_t>(payload[0]));
            const auto count = ReadU32LE(payload.data() + 1);

            TagArray arr{elemType};
            const auto elemSize = arr.ElementSize();
            if (elemSize == 0) {
                throw std::runtime_error{
                    std::format("Tags: unsupported B-array element type '{}'", elemType)};
            }

            const auto dataSize = static_cast<std::size_t>(count) * elemSize;
            const auto expectedSize = static_cast<std::size_t>(5) + dataSize;
            RequirePayloadSize(type, payload, expectedSize);

            arr.Resize(count);
            std::ranges::copy_n(payload.data() + 5, dataSize, arr.MutableData().data());
            return TagValue{std::move(arr)};
        }
        default:
            throw std::runtime_error{std::format("Tags: unsupported BAM tag type '{}'", type)};
    }
}

EncodedTagValuePayload EncodeTagValueToBamPayload(const TagValue& value)
{
    std::vector<std::byte> serialized;
    std::visit(BamSerializeVisitor{serialized}, value);
    if (std::empty(serialized)) {
        throw std::runtime_error{"Tags: empty serialized tag payload"};
    }

    EncodedTagValuePayload result;
    result.Type = static_cast<char>(std::to_integer<std::uint8_t>(serialized.front()));
    serialized.erase(std::begin(serialized));
    result.Payload = std::move(serialized);
    return result;
}

// --- TagKey ---

std::string TagKey::ToString() const { return std::string{First(), Second()}; }

// --- TagArray ---

TagArray::TagArray(char elementType) : elementType_{elementType} {}

char TagArray::ElementType() const { return elementType_; }

std::uint32_t TagArray::Count() const
{
    const std::size_t elemSize{ElementSize()};
    return (elemSize > 0) ? (std::size(data_) / elemSize) : 0;
}

std::size_t TagArray::ElementSize() const
{
    switch (elementType_) {
        case 'c':
        case 'C':
            return 1;
        case 's':
        case 'S':
            return 2;
        case 'i':
        case 'I':
        case 'f':
            return 4;
        default:
            return 0;
    }
}

std::span<const std::byte> TagArray::Data() const { return data_; }

std::span<std::byte> TagArray::MutableData() { return data_; }

void TagArray::Resize(std::uint32_t count) { data_.resize(count * ElementSize()); }

void TagArray::AppendInt8(std::int8_t v) { data_.push_back(static_cast<std::byte>(v)); }

void TagArray::AppendUInt8(std::uint8_t v) { data_.push_back(static_cast<std::byte>(v)); }

void TagArray::AppendInt16(std::int16_t v) { AppendI16LE(data_, v); }

void TagArray::AppendUInt16(std::uint16_t v) { AppendU16LE(data_, v); }

void TagArray::AppendInt32(std::int32_t v) { AppendI32LE(data_, v); }

void TagArray::AppendUInt32(std::uint32_t v) { AppendU32LE(data_, v); }

void TagArray::AppendFloat(float v) { AppendF32LE(data_, v); }

// --- TagMap ---

const TagValue* TagMap::Get(TagKey key) const
{
    const auto it{std::ranges::find_if(entries_, [key](const Entry& e) { return e.first == key; })};
    return (it != std::ranges::end(entries_)) ? &it->second : nullptr;
}

void TagMap::Set(TagKey key, TagValue value)
{
    const auto it{std::ranges::find_if(entries_, [key](const Entry& e) { return e.first == key; })};
    if (it != std::ranges::end(entries_)) {
        it->second = std::move(value);
    } else {
        entries_.emplace_back(key, std::move(value));
    }
}

bool TagMap::Remove(TagKey key)
{
    const auto it{std::ranges::find_if(entries_, [key](const Entry& e) { return e.first == key; })};
    if (it == std::ranges::end(entries_)) {
        return false;
    }
    entries_.erase(it);
    return true;
}

bool TagMap::Contains(TagKey key) const { return Get(key) != nullptr; }

std::span<const TagMap::Entry> TagMap::Entries() const { return entries_; }

std::size_t TagMap::Size() const { return std::size(entries_); }

bool TagMap::Empty() const { return std::empty(entries_); }

void TagMap::Append(TagKey key, TagValue value) { entries_.emplace_back(key, std::move(value)); }

// --- Filters ---

SortedTagKeySet::SortedTagKeySet(std::initializer_list<TagKey> keys) : keys_{keys}
{
    std::ranges::sort(
        keys_, [](const TagKey& lhs, const TagKey& rhs) { return lhs.Value() < rhs.Value(); });
}

bool SortedTagKeySet::Contains(TagKey key) const
{
    const auto it = std::ranges::lower_bound(
        keys_, key, [](const TagKey& lhs, const TagKey& rhs) { return lhs.Value() < rhs.Value(); });
    return (it != std::ranges::end(keys_)) && (*it == key);
}

DropTags::DropTags(std::initializer_list<TagKey> keys) : keys_{keys} {}

bool DropTags::ShouldDrop(TagKey key) const { return keys_.Contains(key); }

KeepTags::KeepTags(std::initializer_list<TagKey> keys) : keys_{keys} {}

bool KeepTags::ShouldKeep(TagKey key) const { return keys_.Contains(key); }

// --- Parsing/serialization ---

TagMap ParseTagsFromBam(std::span<const std::byte> data)
{
    TagMap result;
    std::size_t offset{0};

    while ((offset + 3) <= std::size(data)) {
        const char c1{static_cast<char>(data[offset])};
        const char c2{static_cast<char>(data[offset + 1])};
        const char type{static_cast<char>(data[offset + 2])};
        offset += 3;

        const TagKey key{c1, c2};

        if (type == 'A') {
            if (offset >= std::size(data)) {
                break;
            }
            result.Append(key, TagValue{static_cast<char>(data[offset])});
            offset += 1;
        } else if (type == 'c') {
            if ((offset + 1) > std::size(data)) {
                break;
            }
            const std::int8_t v{static_cast<std::int8_t>(data[offset])};
            result.Append(key, TagValue{std::int64_t{v}});
            offset += 1;
        } else if (type == 'C') {
            if ((offset + 1) > std::size(data)) {
                break;
            }
            const std::uint8_t v{static_cast<std::uint8_t>(data[offset])};
            result.Append(key, TagValue{std::int64_t{v}});
            offset += 1;
        } else if (type == 's') {
            if ((offset + 2) > std::size(data)) {
                break;
            }
            const std::int16_t v{ReadI16LE(std::data(data) + offset)};
            result.Append(key, TagValue{std::int64_t{v}});
            offset += 2;
        } else if (type == 'S') {
            if ((offset + 2) > std::size(data)) {
                break;
            }
            const std::uint16_t v{ReadU16LE(std::data(data) + offset)};
            result.Append(key, TagValue{std::int64_t{v}});
            offset += 2;
        } else if (type == 'i') {
            if ((offset + 4) > std::size(data)) {
                break;
            }
            const std::int32_t v{ReadI32LE(std::data(data) + offset)};
            result.Append(key, TagValue{std::int64_t{v}});
            offset += 4;
        } else if (type == 'I') {
            if ((offset + 4) > std::size(data)) {
                break;
            }
            const std::uint32_t v{ReadU32LE(std::data(data) + offset)};
            result.Append(key, TagValue{std::int64_t{v}});
            offset += 4;
        } else if (type == 'f') {
            if ((offset + 4) > std::size(data)) {
                break;
            }
            const float v{ReadF32LE(std::data(data) + offset)};
            result.Append(key, TagValue{v});
            offset += 4;
        } else if (type == 'Z') {
            const std::size_t start{offset};
            while (offset < std::size(data) && data[offset] != std::byte{0}) {
                ++offset;
            }
            std::string s(reinterpret_cast<const char*>(std::data(data) + start), offset - start);
            result.Append(key, TagValue{std::move(s)});
            if (offset < std::size(data)) {
                ++offset;  // skip NUL
            }
        } else if (type == 'H') {
            const std::size_t start{offset};
            while (offset < std::size(data) && data[offset] != std::byte{0}) {
                ++offset;
            }
            std::string s(reinterpret_cast<const char*>(std::data(data) + start), offset - start);
            result.Append(key, TagValue{HexString{std::move(s)}});
            if (offset < std::size(data)) {
                ++offset;
            }
        } else if (type == 'B') {
            if ((offset + 5) > std::size(data)) {
                break;
            }
            const char elemType{static_cast<char>(data[offset])};
            offset += 1;
            const std::uint32_t count{ReadU32LE(std::data(data) + offset)};
            offset += 4;

            TagArray arr{elemType};
            const std::size_t elemSize{arr.ElementSize()};
            const std::size_t totalBytes{count * elemSize};
            if ((offset + totalBytes) > std::size(data)) {
                break;
            }
            arr.Resize(count);
            std::ranges::copy_n(std::data(data) + offset, totalBytes, std::data(arr.MutableData()));
            offset += totalBytes;
            result.Append(key, TagValue{std::move(arr)});
        }
    }

    return result;
}

std::optional<std::pair<TagKey, TagValue>> ParseTagFromSam(std::string_view text)
{
    // Format: XX:T:VALUE  (minimum length 5: 2 key + ':' + type + ':')
    if (std::size(text) < 5) {
        return std::nullopt;
    }
    if ((text[2] != ':') || (text[4] != ':')) {
        return std::nullopt;
    }

    const TagKey key{text[0], text[1]};
    const char type{text[3]};
    const std::string_view valueStr{text.substr(5)};

    if (type == 'A') {
        if (std::size(valueStr) != 1) {
            return std::nullopt;
        }
        return std::pair{key, TagValue{valueStr[0]}};
    } else if (type == 'i') {
        std::int64_t v{0};
        const std::from_chars_result parseResult{
            std::from_chars(std::data(valueStr), std::data(valueStr) + std::size(valueStr), v)};
        if ((parseResult.ec != std::errc{}) ||
            (parseResult.ptr != (std::data(valueStr) + std::size(valueStr)))) {
            return std::nullopt;
        }
        return std::pair{key, TagValue{v}};
    } else if (type == 'f') {
        float v{0.0f};
        const std::from_chars_result parseResult{
            std::from_chars(std::data(valueStr), std::data(valueStr) + std::size(valueStr), v)};
        if ((parseResult.ec != std::errc{}) ||
            (parseResult.ptr != (std::data(valueStr) + std::size(valueStr)))) {
            return std::nullopt;
        }
        return std::pair{key, TagValue{v}};
    } else if (type == 'Z') {
        return std::pair{key, TagValue{std::string{valueStr}}};
    } else if (type == 'H') {
        if ((std::size(valueStr) % 2) != 0) {
            return std::nullopt;
        }
        if (!std::ranges::all_of(valueStr, IsHexDigit)) {
            return std::nullopt;
        }
        return std::pair{key, TagValue{HexString{std::string{valueStr}}}};
    } else if (type == 'B') {
        // Format: B:T,v1,v2,...
        if (std::size(valueStr) < 3) {
            return std::nullopt;
        }
        const char elemType{valueStr[0]};
        if ((elemType != 'c') && (elemType != 'C') && (elemType != 's') && (elemType != 'S') &&
            (elemType != 'i') && (elemType != 'I') && (elemType != 'f')) {
            return std::nullopt;
        }
        if (valueStr[1] != ',') {
            return std::nullopt;
        }

        TagArray arr{elemType};

        // Parse comma-separated values after element type
        std::size_t pos{2};  // skip element type + comma
        while (pos < std::size(valueStr)) {
            // Find end of this value
            const std::size_t commaPos{valueStr.find(',', pos)};
            const std::size_t endPos{(commaPos != std::string_view::npos) ? commaPos
                                                                          : std::size(valueStr)};
            const std::string_view elem{valueStr.substr(pos, endPos - pos)};
            if (std::empty(elem)) {
                return std::nullopt;
            }

            if (elemType == 'f') {
                float v{0.0f};
                const std::from_chars_result parseResult{
                    std::from_chars(std::data(elem), std::data(elem) + std::size(elem), v)};
                if ((parseResult.ec != std::errc{}) ||
                    (parseResult.ptr != (std::data(elem) + std::size(elem)))) {
                    return std::nullopt;
                }
                arr.AppendFloat(v);
            } else {
                // Integer element types
                std::int64_t v{0};
                const std::from_chars_result parseResult{
                    std::from_chars(std::data(elem), std::data(elem) + std::size(elem), v)};
                if ((parseResult.ec != std::errc{}) ||
                    (parseResult.ptr != (std::data(elem) + std::size(elem)))) {
                    return std::nullopt;
                }

                switch (elemType) {
                    case 'c':
                        if ((v < std::numeric_limits<std::int8_t>::min()) ||
                            (v > std::numeric_limits<std::int8_t>::max())) {
                            return std::nullopt;
                        }
                        arr.AppendInt8(v);
                        break;
                    case 'C':
                        if ((v < std::numeric_limits<std::uint8_t>::min()) ||
                            (v > std::numeric_limits<std::uint8_t>::max())) {
                            return std::nullopt;
                        }
                        arr.AppendUInt8(v);
                        break;
                    case 's':
                        if ((v < std::numeric_limits<std::int16_t>::min()) ||
                            (v > std::numeric_limits<std::int16_t>::max())) {
                            return std::nullopt;
                        }
                        arr.AppendInt16(v);
                        break;
                    case 'S':
                        if ((v < std::numeric_limits<std::uint16_t>::min()) ||
                            (v > std::numeric_limits<std::uint16_t>::max())) {
                            return std::nullopt;
                        }
                        arr.AppendUInt16(v);
                        break;
                    case 'i':
                        if ((v < std::numeric_limits<std::int32_t>::min()) ||
                            (v > std::numeric_limits<std::int32_t>::max())) {
                            return std::nullopt;
                        }
                        arr.AppendInt32(v);
                        break;
                    case 'I':
                        if ((v < std::numeric_limits<std::uint32_t>::min()) ||
                            (v > std::numeric_limits<std::uint32_t>::max())) {
                            return std::nullopt;
                        }
                        arr.AppendUInt32(v);
                        break;
                    default:
                        return std::nullopt;
                }
            }
            pos = (commaPos == std::string_view::npos) ? endPos : (endPos + 1);
        }

        return std::pair{key, TagValue{std::move(arr)}};
    }

    return std::nullopt;
}

std::vector<std::byte> SerializeTagsToBam(const TagMap& tags)
{
    std::vector<std::byte> result;

    for (const auto& [key, value] : tags.Entries()) {
        // Write 2-byte key
        result.push_back(static_cast<std::byte>(key.First()));
        result.push_back(static_cast<std::byte>(key.Second()));

        std::visit(BamSerializeVisitor{result}, value);
    }

    return result;
}

std::string SerializeTagToSam(TagKey key, const TagValue& value)
{
    std::string result;
    result += key.First();
    result += key.Second();
    result += ':';

    std::visit(SamSerializeVisitor{result}, value);

    return result;
}

namespace {

void AppendIntRaw(std::string& out, std::int64_t v)
{
    std::array<char, 24> buf{};
    const auto [ptr, ec]{std::to_chars(std::data(buf), std::data(buf) + std::size(buf), v)};
    out.append(std::data(buf), ptr);
}

void AppendFloatRaw(std::string& out, float v)
{
    std::array<char, 32> buf{};
    const int n{std::snprintf(std::data(buf), std::size(buf), "%g", static_cast<double>(v))};
    if (n > 0) {
        const std::size_t len = n;
        out.append(std::data(buf), len);
    }
}

/// \brief Fast-path: serialize B:C array using resize_and_overwrite (no zero-fill).
void SerializeBArrayUInt8(const std::byte* data, std::uint32_t count, std::string& out)
{
    const std::size_t startPos{std::size(out)};
    const std::size_t maxChars{count * 4U};

    out.resize_and_overwrite(
        startPos + maxChars,
        [data, count, startPos](char* buf, std::size_t /*bufSize*/) -> std::size_t {
            char* dest{buf + startPos};
            const auto appendChar = [&dest](const char c) {
                *dest = c;
                ++dest;
            };
            for (std::uint32_t i{0}; i < count; ++i) {
                const unsigned v = static_cast<std::uint8_t>(data[i]);
                appendChar(',');
                if (v >= 100) {
                    appendChar('0' + v / 100);
                    appendChar('0' + (v / 10) % 10);
                    appendChar('0' + v % 10);
                } else if (v >= 10) {
                    appendChar('0' + v / 10);
                    appendChar('0' + v % 10);
                } else {
                    appendChar('0' + v);
                }
            }
            return dest - buf;
        });
}

/// \brief Fast-path: serialize B:c array using resize_and_overwrite (no zero-fill).
void SerializeBArrayInt8(const std::byte* data, std::uint32_t count, std::string& out)
{
    const std::size_t startPos{std::size(out)};
    const std::size_t maxChars{count * 5U};

    out.resize_and_overwrite(
        startPos + maxChars,
        [data, count, startPos](char* buf, std::size_t /*bufSize*/) -> std::size_t {
            char* dest{buf + startPos};
            const auto appendChar = [&dest](const char c) {
                *dest = c;
                ++dest;
            };
            for (std::uint32_t i{0}; i < count; ++i) {
                const int v{static_cast<std::int8_t>(data[i])};
                appendChar(',');
                const int absV{(v < 0) ? -v : v};
                if (v < 0) {
                    appendChar('-');
                }
                if (absV >= 100) {
                    appendChar('0' + absV / 100);
                    appendChar('0' + (absV / 10) % 10);
                    appendChar('0' + absV % 10);
                } else if (absV >= 10) {
                    appendChar('0' + absV / 10);
                    appendChar('0' + absV % 10);
                } else {
                    appendChar('0' + absV);
                }
            }
            return dest - buf;
        });
}

/// \brief Serialize B:S/s/I/i array elements using to_chars with correct type width.
template <typename T>
void SerializeBArrayInt(const std::byte* data, std::uint32_t count, std::string& out)
{
    std::array<char, 16> buf{};
    for (std::uint32_t i{0}; i < count; ++i) {
        out += ',';
        T v{};
        if constexpr (std::is_same_v<T, std::int16_t>) {
            v = ReadI16LE(data + i * sizeof(T));
        } else if constexpr (std::is_same_v<T, std::uint16_t>) {
            v = ReadU16LE(data + i * sizeof(T));
        } else if constexpr (std::is_same_v<T, std::int32_t>) {
            v = ReadI32LE(data + i * sizeof(T));
        } else if constexpr (std::is_same_v<T, std::uint32_t>) {
            v = ReadU32LE(data + i * sizeof(T));
        }
        const auto [ptr, ec]{std::to_chars(std::data(buf), std::data(buf) + std::size(buf), v)};
        out.append(std::data(buf), ptr);
    }
}

void SerializeBArrayFloat(const std::byte* data, std::uint32_t count, std::string& out)
{
    for (std::uint32_t i{0}; i < count; ++i) {
        out += ',';
        AppendFloatRaw(out, ReadF32LE(data + i * 4));
    }
}

}  // namespace

void SerializeRawTagsToSam(std::span<const std::byte> data, std::string& out)
{
    std::size_t offset{0};

    while ((offset + 3) <= std::size(data)) {
        const char c1{static_cast<char>(data[offset])};
        const char c2{static_cast<char>(data[offset + 1])};
        const char type{static_cast<char>(data[offset + 2])};
        offset += 3;

        out += '\t';
        out += c1;
        out += c2;
        out += ':';

        if (type == 'A') {
            if (offset >= std::size(data)) {
                break;
            }
            out += "A:";
            out += static_cast<char>(data[offset]);
            offset += 1;
        } else if (type == 'c') {
            if ((offset + 1) > std::size(data)) {
                break;
            }
            out += "i:";
            AppendIntRaw(out, static_cast<std::int8_t>(data[offset]));
            offset += 1;
        } else if (type == 'C') {
            if ((offset + 1) > std::size(data)) {
                break;
            }
            out += "i:";
            AppendIntRaw(out, static_cast<std::uint8_t>(data[offset]));
            offset += 1;
        } else if (type == 's') {
            if ((offset + 2) > std::size(data)) {
                break;
            }
            out += "i:";
            AppendIntRaw(out, ReadI16LE(std::data(data) + offset));
            offset += 2;
        } else if (type == 'S') {
            if ((offset + 2) > std::size(data)) {
                break;
            }
            out += "i:";
            AppendIntRaw(out, ReadU16LE(std::data(data) + offset));
            offset += 2;
        } else if (type == 'i') {
            if ((offset + 4) > std::size(data)) {
                break;
            }
            out += "i:";
            AppendIntRaw(out, ReadI32LE(std::data(data) + offset));
            offset += 4;
        } else if (type == 'I') {
            if ((offset + 4) > std::size(data)) {
                break;
            }
            out += "i:";
            AppendIntRaw(out, ReadU32LE(std::data(data) + offset));
            offset += 4;
        } else if (type == 'f') {
            if ((offset + 4) > std::size(data)) {
                break;
            }
            out += "f:";
            AppendFloatRaw(out, ReadF32LE(std::data(data) + offset));
            offset += 4;
        } else if (type == 'Z') {
            out += "Z:";
            while ((offset < std::size(data)) && (data[offset] != std::byte{0})) {
                out += static_cast<char>(data[offset]);
                ++offset;
            }
            if (offset < std::size(data)) {
                ++offset;  // skip NUL
            }
        } else if (type == 'H') {
            out += "H:";
            while ((offset < std::size(data)) && (data[offset] != std::byte{0})) {
                out += static_cast<char>(data[offset]);
                ++offset;
            }
            if (offset < std::size(data)) {
                ++offset;
            }
        } else if (type == 'B') {
            if ((offset + 5) > std::size(data)) {
                break;
            }
            const char elemType{static_cast<char>(data[offset])};
            offset += 1;
            const std::uint32_t count{ReadU32LE(std::data(data) + offset)};
            offset += 4;

            out += "B:";
            out += elemType;

            // Fast-path for byte arrays (dominant in PacBio kinetics data)
            if (elemType == 'C') {
                if ((offset + count) > std::size(data)) {
                    break;
                }
                SerializeBArrayUInt8(std::data(data) + offset, count, out);
                offset += count;
            } else if (elemType == 'c') {
                if ((offset + count) > std::size(data)) {
                    break;
                }
                SerializeBArrayInt8(std::data(data) + offset, count, out);
                offset += count;
            } else if (elemType == 's') {
                if ((offset + count * 2) > std::size(data)) {
                    break;
                }
                SerializeBArrayInt<std::int16_t>(std::data(data) + offset, count, out);
                offset += count * 2;
            } else if (elemType == 'S') {
                if ((offset + count * 2) > std::size(data)) {
                    break;
                }
                SerializeBArrayInt<std::uint16_t>(std::data(data) + offset, count, out);
                offset += count * 2;
            } else if (elemType == 'i') {
                if ((offset + count * 4) > std::size(data)) {
                    break;
                }
                SerializeBArrayInt<std::int32_t>(std::data(data) + offset, count, out);
                offset += count * 4;
            } else if (elemType == 'I') {
                if ((offset + count * 4) > std::size(data)) {
                    break;
                }
                SerializeBArrayInt<std::uint32_t>(std::data(data) + offset, count, out);
                offset += count * 4;
            } else if (elemType == 'f') {
                if ((offset + count * 4) > std::size(data)) {
                    break;
                }
                SerializeBArrayFloat(std::data(data) + offset, count, out);
                offset += count * 4;
            }
        }
    }
}

// --- Zero-allocation tag serialization ---

namespace {

struct BamSizeVisitor
{
    std::size_t operator()(char /*v*/) const
    {
        // type byte + value byte
        return 2;
    }

    std::size_t operator()(std::int64_t v) const
    {
        const char bamType{SmallestIntType(v)};
        // type byte + value bytes
        switch (bamType) {
            case 'c':
            case 'C':
                return 2;
            case 's':
            case 'S':
                return 3;
            case 'i':
            case 'I':
                return 5;
            default:
                return 1;
        }
    }

    std::size_t operator()(float /*v*/) const
    {
        // type byte + 4 value bytes
        return 5;
    }

    std::size_t operator()(std::string_view v) const
    {
        // type byte + string chars + NUL
        return 1 + std::size(v) + 1;
    }

    std::size_t operator()(const HexString& v) const
    {
        // type byte + hex chars + NUL
        return 1 + std::size(v.value) + 1;
    }

    std::size_t operator()(const TagArray& v) const
    {
        // type byte + element type byte + 4-byte count + element data
        return 1 + 1 + 4 + std::size(v.Data());
    }
};

void WriteU16Raw(std::byte*& dest, std::uint16_t v)
{
    *dest++ = static_cast<std::byte>(v & 0xFFU);
    *dest++ = static_cast<std::byte>((v >> 8U) & 0xFFU);
}

void WriteU32Raw(std::byte*& dest, std::uint32_t v)
{
    *dest++ = static_cast<std::byte>(v & 0xFFU);
    *dest++ = static_cast<std::byte>((v >> 8U) & 0xFFU);
    *dest++ = static_cast<std::byte>((v >> 16U) & 0xFFU);
    *dest++ = static_cast<std::byte>((v >> 24U) & 0xFFU);
}

void WriteI16Raw(std::byte*& dest, std::int16_t v)
{
    WriteU16Raw(dest, std::bit_cast<std::uint16_t>(v));
}

void WriteI32Raw(std::byte*& dest, std::int32_t v)
{
    WriteU32Raw(dest, std::bit_cast<std::uint32_t>(v));
}

void WriteF32Raw(std::byte*& dest, float v) { WriteU32Raw(dest, std::bit_cast<std::uint32_t>(v)); }

struct BamAppendVisitor
{
    std::byte*& dest;

    void AppendByte(const std::byte v) const
    {
        *dest = v;
        ++dest;
    }

    void operator()(char v) const
    {
        AppendByte(static_cast<std::byte>('A'));
        AppendByte(static_cast<std::byte>(v));
    }

    void operator()(std::int64_t v) const
    {
        const char bamType{SmallestIntType(v)};
        AppendByte(static_cast<std::byte>(bamType));
        switch (bamType) {
            case 'c': {
                const std::int8_t sv = v;
                AppendByte(static_cast<std::byte>(sv));
                break;
            }
            case 'C': {
                const std::uint8_t sv = v;
                AppendByte(static_cast<std::byte>(sv));
                break;
            }
            case 's': {
                const std::int16_t sv = v;
                WriteI16Raw(dest, sv);
                break;
            }
            case 'S': {
                const std::uint16_t sv = v;
                WriteU16Raw(dest, sv);
                break;
            }
            case 'i': {
                const std::int32_t sv = v;
                WriteI32Raw(dest, sv);
                break;
            }
            case 'I': {
                const std::uint32_t sv = v;
                WriteU32Raw(dest, sv);
                break;
            }
            default:
                break;
        }
    }

    void operator()(float v) const
    {
        *dest++ = static_cast<std::byte>('f');
        WriteF32Raw(dest, v);
    }

    void operator()(std::string_view v) const
    {
        AppendByte(static_cast<std::byte>('Z'));
        for (const char ch : v) {
            AppendByte(static_cast<std::byte>(ch));
        }
        AppendByte(std::byte{0});
    }

    void operator()(const HexString& v) const
    {
        AppendByte(static_cast<std::byte>('H'));
        for (const char ch : v.value) {
            AppendByte(static_cast<std::byte>(ch));
        }
        AppendByte(std::byte{0});
    }

    void operator()(const TagArray& v) const
    {
        AppendByte(static_cast<std::byte>('B'));
        AppendByte(static_cast<std::byte>(v.ElementType()));
        const std::uint32_t count{v.Count()};
        WriteU32Raw(dest, count);
        const std::span<const std::byte> data{v.Data()};
        std::ranges::copy_n(std::data(data), std::size(data), dest);
        dest += std::size(data);
    }
};

}  // namespace

std::size_t SerializedBamSize(const TagMap& tags)
{
    std::size_t total{0};
    for (const auto& [key, value] : tags.Entries()) {
        // 2-byte key + value size (includes type byte)
        total += 2 + std::visit(BamSizeVisitor{}, value);
    }
    return total;
}

void AppendTagsToBam(const TagMap& tags, std::byte* dest)
{
    const auto appendByte = [&dest](const std::byte value) {
        *dest = value;
        ++dest;
    };
    for (const auto& [key, value] : tags.Entries()) {
        appendByte(static_cast<std::byte>(key.First()));
        appendByte(static_cast<std::byte>(key.Second()));
        std::visit(BamAppendVisitor{dest}, value);
    }
}

}  // namespace Samoa
}  // namespace PacBio
