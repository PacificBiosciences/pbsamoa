#include <pbsamoa/core/Tags.hpp>

#include "CramInternal.hpp"

#include <pbsamoa/core/Endian.hpp>

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

std::size_t TagArrayElementSize(char elementType)
{
    switch (elementType) {
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

void WriteU16Raw(std::byte*& dest, std::uint16_t v);
void WriteU32Raw(std::byte*& dest, std::uint32_t v);
void WriteI16Raw(std::byte*& dest, std::int16_t v);
void WriteI32Raw(std::byte*& dest, std::int32_t v);

void AppendIntegerTagValue(std::vector<std::byte>& result, char bamType, std::int64_t value)
{
    switch (bamType) {
        case 'c':
        case 'C':
            result.push_back(static_cast<std::byte>(value));
            break;
        case 's': {
            const auto narrowed = static_cast<std::int16_t>(value);
            AppendI16LE(result, narrowed);
            break;
        }
        case 'S': {
            const auto narrowed = static_cast<std::uint16_t>(value);
            AppendU16LE(result, narrowed);
            break;
        }
        case 'i': {
            const auto narrowed = static_cast<std::int32_t>(value);
            AppendI32LE(result, narrowed);
            break;
        }
        case 'I': {
            const auto narrowed = static_cast<std::uint32_t>(value);
            AppendU32LE(result, narrowed);
            break;
        }
        default:
            break;
    }
}

void AppendRawByte(std::byte*& dest, std::byte value)
{
    *dest = value;
    ++dest;
}

void AppendIntegerTagValue(std::byte*& dest, char bamType, std::int64_t value)
{
    switch (bamType) {
        case 'c': {
            const auto narrowed = static_cast<std::int8_t>(value);
            AppendRawByte(dest, static_cast<std::byte>(narrowed));
            break;
        }
        case 'C': {
            const auto narrowed = static_cast<std::uint8_t>(value);
            AppendRawByte(dest, static_cast<std::byte>(narrowed));
            break;
        }
        case 's': {
            const auto narrowed = static_cast<std::int16_t>(value);
            WriteI16Raw(dest, narrowed);
            break;
        }
        case 'S': {
            const auto narrowed = static_cast<std::uint16_t>(value);
            WriteU16Raw(dest, narrowed);
            break;
        }
        case 'i': {
            const auto narrowed = static_cast<std::int32_t>(value);
            WriteI32Raw(dest, narrowed);
            break;
        }
        case 'I': {
            const auto narrowed = static_cast<std::uint32_t>(value);
            WriteU32Raw(dest, narrowed);
            break;
        }
        default:
            break;
    }
}

bool AppendIntegerArrayElement(TagArray& array, char elemType, std::int64_t value)
{
    switch (elemType) {
        case 'c':
            if ((value < std::numeric_limits<std::int8_t>::min()) ||
                (value > std::numeric_limits<std::int8_t>::max())) {
                return false;
            }
            array.AppendInt8(value);
            return true;
        case 'C':
            if ((value < std::numeric_limits<std::uint8_t>::min()) ||
                (value > std::numeric_limits<std::uint8_t>::max())) {
                return false;
            }
            array.AppendUInt8(value);
            return true;
        case 's':
            if ((value < std::numeric_limits<std::int16_t>::min()) ||
                (value > std::numeric_limits<std::int16_t>::max())) {
                return false;
            }
            array.AppendInt16(value);
            return true;
        case 'S':
            if ((value < std::numeric_limits<std::uint16_t>::min()) ||
                (value > std::numeric_limits<std::uint16_t>::max())) {
                return false;
            }
            array.AppendUInt16(value);
            return true;
        case 'i':
            if ((value < std::numeric_limits<std::int32_t>::min()) ||
                (value > std::numeric_limits<std::int32_t>::max())) {
                return false;
            }
            array.AppendInt32(value);
            return true;
        case 'I':
            if ((value < std::numeric_limits<std::uint32_t>::min()) ||
                (value > std::numeric_limits<std::uint32_t>::max())) {
                return false;
            }
            array.AppendUInt32(value);
            return true;
        default:
            return false;
    }
}

void AppendTagArrayElementToSam(std::string& result, char elemType, const std::byte* data)
{
    switch (elemType) {
        case 'f':
            std::format_to(std::back_inserter(result), "{:g}", ReadF32LE(data));
            break;
        case 'c':
            std::format_to(std::back_inserter(result), "{}",
                           static_cast<std::int8_t>(std::to_integer<std::uint8_t>(*data)));
            break;
        case 'C':
            std::format_to(std::back_inserter(result), "{}", std::to_integer<std::uint8_t>(*data));
            break;
        case 's':
            std::format_to(std::back_inserter(result), "{}", ReadI16LE(data));
            break;
        case 'S':
            std::format_to(std::back_inserter(result), "{}", ReadU16LE(data));
            break;
        case 'i':
            std::format_to(std::back_inserter(result), "{}", ReadI32LE(data));
            break;
        case 'I':
            std::format_to(std::back_inserter(result), "{}", ReadU32LE(data));
            break;
        default:
            break;
    }
}

char SmallestIntType(std::int64_t v)
{
    // A BAM integer tag is representable only in [INT32_MIN, UINT32_MAX]; htslib
    // (bam_aux_update_int) rejects wider values. Fail loud instead of truncating to 32 bits.
    if ((v < std::numeric_limits<std::int32_t>::min()) ||
        (v > std::numeric_limits<std::uint32_t>::max())) {
        throw std::runtime_error{
            std::format("Tags: integer tag value {} out of BAM-representable 32-bit range", v)};
    }
    if ((v >= 0) && (v <= std::numeric_limits<std::uint8_t>::max())) {
        return 'C';
    }
    if ((v >= std::numeric_limits<std::int8_t>::min()) &&
        (v <= std::numeric_limits<std::int8_t>::max())) {
        return 'c';
    }
    if ((v >= 0) && (v <= std::numeric_limits<std::uint16_t>::max())) {
        return 'S';
    }
    if ((v >= std::numeric_limits<std::int16_t>::min()) &&
        (v <= std::numeric_limits<std::int16_t>::max())) {
        return 's';
    }
    if ((v >= 0) && (v <= std::numeric_limits<std::uint32_t>::max())) {
        return 'I';
    }
    return 'i';
}

bool IsHexDigit(char c)
{
    return ((c >= '0') && (c <= '9')) || ((c >= 'a') && (c <= 'f')) || ((c >= 'A') && (c <= 'F'));
}

template <typename T>
bool ParseExact(std::string_view text, T& value)
{
    const char* const begin{std::data(text)};
    const char* const end{begin + std::size(text)};
    const auto [ptr, ec]{std::from_chars(begin, end, value)};
    return (ec == std::errc{}) && (ptr == end);
}

std::string ReadTextPayload(std::span<const std::byte> payload)
{
    if (!std::empty(payload) && (payload.back() == std::byte{0})) {
        payload = payload.first(std::size(payload) - 1);
    }
    return std::string{reinterpret_cast<const char*>(payload.data()), std::size(payload)};
}

std::string ReadNulTerminatedText(std::span<const std::byte> data, std::size_t& offset)
{
    const std::size_t start{offset};
    while ((offset < std::size(data)) && (data[offset] != std::byte{0})) {
        ++offset;
    }

    if (offset == std::size(data)) {
        throw std::runtime_error{"Tags: unterminated BAM auxiliary text payload"};
    }
    std::string text{reinterpret_cast<const char*>(std::data(data) + start), offset - start};
    ++offset;
    return text;
}

void RequireBamAuxBytes(std::span<const std::byte> data, std::size_t offset, std::size_t count)
{
    if ((offset > std::size(data)) || (count > (std::size(data) - offset))) {
        throw std::runtime_error{"Tags: truncated BAM auxiliary data"};
    }
}

void RequirePayloadSize(char type, std::span<const std::byte> payload, std::size_t expected)
{
    if (std::size(payload) != expected) {
        throw std::runtime_error{std::format("Tags: type '{}' requires {} payload byte(s), got {}",
                                             type, expected, std::size(payload))};
    }
}

using EntryIter = std::vector<TagMap::Entry>::iterator;
using ConstEntryIter = std::vector<TagMap::Entry>::const_iterator;

EntryIter LowerBoundEntry(std::vector<TagMap::Entry>& entries, TagKey key)
{
    return std::ranges::lower_bound(entries, key.Value(), {},
                                    [](const TagMap::Entry& entry) { return entry.first.Value(); });
}

ConstEntryIter LowerBoundEntry(const std::vector<TagMap::Entry>& entries, TagKey key)
{
    return std::ranges::lower_bound(entries, key.Value(), {},
                                    [](const TagMap::Entry& entry) { return entry.first.Value(); });
}

EntryIter FindEntry(std::vector<TagMap::Entry>& entries, TagKey key)
{
    const EntryIter it{LowerBoundEntry(entries, key)};
    if ((it == std::end(entries)) || (it->first != key)) {
        return std::end(entries);
    }
    return it;
}

ConstEntryIter FindEntry(const std::vector<TagMap::Entry>& entries, TagKey key)
{
    const ConstEntryIter it{LowerBoundEntry(entries, key)};
    if ((it == std::cend(entries)) || (it->first != key)) {
        return std::cend(entries);
    }
    return it;
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
        AppendIntegerTagValue(result, bamType, v);
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
        const std::size_t elemSize{TagArrayElementSize(elemType)};
        const std::byte* element{std::data(data)};

        for (std::size_t i{0}; i < count; ++i) {
            result += ',';
            AppendTagArrayElementToSam(result, elemType, element);
            element += elemSize;
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
            return TagValue{ReadTextPayload(payload)};
        }
        case 'H': {
            // SAMv1 §4.2.4: 'H' is a byte array encoded as an even-length run of hex
            // digits. Validate here so the BAM decode path matches the SAM-text path and
            // never round-trips malformed hex back out as non-conformant SAM.
            std::string hex{ReadTextPayload(payload)};
            if (((std::size(hex) % 2) != 0) || !std::ranges::all_of(hex, IsHexDigit)) {
                throw std::runtime_error{"Tags: type 'H' payload is not valid hex"};
            }
            return TagValue{HexString{std::move(hex)}};
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
    if (elemSize == 0) {
        return 0;
    }
    return std::size(data_) / elemSize;
}

std::size_t TagArray::ElementSize() const { return TagArrayElementSize(elementType_); }

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
    const ConstEntryIter it{FindEntry(entries_, key)};
    if (it == std::cend(entries_)) {
        return nullptr;
    }
    return &it->second;
}

void TagMap::Set(TagKey key, TagValue value)
{
    const EntryIter it{LowerBoundEntry(entries_, key)};
    if ((it != std::end(entries_)) && (it->first == key)) {
        it->second = std::move(value);
        return;
    }
    entries_.insert(it, Entry{key, std::move(value)});
}

bool TagMap::Remove(TagKey key)
{
    const EntryIter it{FindEntry(entries_, key)};
    if (it == std::end(entries_)) {
        return false;
    }
    entries_.erase(it);
    return true;
}

bool TagMap::Contains(TagKey key) const { return FindEntry(entries_, key) != std::cend(entries_); }

std::span<const TagMap::Entry> TagMap::Entries() const { return entries_; }

std::size_t TagMap::Size() const { return std::size(entries_); }

bool TagMap::Empty() const { return std::empty(entries_); }

void TagMap::Append(TagKey key, TagValue value)
{
    entries_.insert(LowerBoundEntry(entries_, key), Entry{key, std::move(value)});
}

// --- Filters ---

SortedTagKeySet::SortedTagKeySet(std::initializer_list<TagKey> keys) : keys_{keys}
{
    std::ranges::sort(keys_, {}, [](TagKey k) { return k.Value(); });
}

bool SortedTagKeySet::Contains(TagKey key) const
{
    const std::uint16_t keyValue{key.Value()};
    const auto it =
        std::ranges::lower_bound(keys_, keyValue, {}, [](TagKey k) { return k.Value(); });
    return (it != std::ranges::end(keys_)) && (it->Value() == keyValue);
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

    while (offset < std::size(data)) {
        RequireBamAuxBytes(data, offset, 3);
        const char c1{static_cast<char>(data[offset])};
        const char c2{static_cast<char>(data[offset + 1])};
        const char type{static_cast<char>(data[offset + 2])};
        offset += 3;

        const TagKey key{c1, c2};

        switch (type) {
            case 'A':
                RequireBamAuxBytes(data, offset, 1);
                result.Append(key, TagValue{static_cast<char>(data[offset])});
                ++offset;
                break;
            case 'c':
                RequireBamAuxBytes(data, offset, 1);
                result.Append(key, TagValue{std::int64_t{static_cast<std::int8_t>(data[offset])}});
                ++offset;
                break;
            case 'C':
                RequireBamAuxBytes(data, offset, 1);
                result.Append(key, TagValue{std::int64_t{static_cast<std::uint8_t>(data[offset])}});
                ++offset;
                break;
            case 's':
                RequireBamAuxBytes(data, offset, 2);
                result.Append(key, TagValue{std::int64_t{ReadI16LE(std::data(data) + offset)}});
                offset += 2;
                break;
            case 'S':
                RequireBamAuxBytes(data, offset, 2);
                result.Append(key, TagValue{std::int64_t{ReadU16LE(std::data(data) + offset)}});
                offset += 2;
                break;
            case 'i':
                RequireBamAuxBytes(data, offset, 4);
                result.Append(key, TagValue{std::int64_t{ReadI32LE(std::data(data) + offset)}});
                offset += 4;
                break;
            case 'I':
                RequireBamAuxBytes(data, offset, 4);
                result.Append(key, TagValue{std::int64_t{ReadU32LE(std::data(data) + offset)}});
                offset += 4;
                break;
            case 'f':
                RequireBamAuxBytes(data, offset, 4);
                result.Append(key, TagValue{ReadF32LE(std::data(data) + offset)});
                offset += 4;
                break;
            case 'Z': {
                result.Append(key, TagValue{ReadNulTerminatedText(data, offset)});
                break;
            }
            case 'H': {
                // SAMv1 §4.2.4: 'H' must be an even-length run of hex digits (matches the
                // SAM-text and single-value decode paths; reject malformed hex loudly).
                std::string hex{ReadNulTerminatedText(data, offset)};
                if (((std::size(hex) % 2) != 0) || !std::ranges::all_of(hex, IsHexDigit)) {
                    throw std::runtime_error{"Tags: type 'H' payload is not valid hex"};
                }
                result.Append(key, TagValue{HexString{std::move(hex)}});
                break;
            }
            case 'B': {
                RequireBamAuxBytes(data, offset, 5);
                const char elemType{static_cast<char>(data[offset])};
                offset += 1;
                const std::uint32_t count{ReadU32LE(std::data(data) + offset)};
                offset += 4;

                TagArray arr{elemType};
                const std::size_t elemSize{arr.ElementSize()};
                if (elemSize == 0) {
                    throw std::runtime_error{
                        std::format("Tags: unsupported B-array element type '{}'", elemType)};
                }
                const std::size_t totalBytes{count * elemSize};
                RequireBamAuxBytes(data, offset, totalBytes);
                arr.Resize(count);
                std::ranges::copy_n(std::data(data) + offset, totalBytes,
                                    std::data(arr.MutableData()));
                offset += totalBytes;
                result.Append(key, TagValue{std::move(arr)});
                break;
            }
            default:
                throw std::runtime_error{std::format("Tags: unsupported BAM tag type '{}'", type)};
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

    // SAMv1 §1.5: TAG must match [A-Za-z][A-Za-z0-9]
    {
        const auto c0{static_cast<unsigned char>(text[0])};
        const auto c1{static_cast<unsigned char>(text[1])};
        const bool firstOk{(c0 >= 'A' && c0 <= 'Z') || (c0 >= 'a' && c0 <= 'z')};
        const bool secondOk{(c1 >= 'A' && c1 <= 'Z') || (c1 >= 'a' && c1 <= 'z') ||
                            (c1 >= '0' && c1 <= '9')};
        if (!firstOk || !secondOk) {
            return std::nullopt;
        }
    }
    const TagKey key{text[0], text[1]};
    const char type{text[3]};
    const std::string_view valueStr{text.substr(5)};

    switch (type) {
        case 'A':
            if (std::size(valueStr) != 1) {
                return std::nullopt;
            }
            // SAMv1 §1.5 type 'A': value must be a printable character [!-~]
            {
                const auto uch{static_cast<unsigned char>(valueStr[0])};
                if (uch < 0x21U || uch > 0x7EU) {
                    return std::nullopt;
                }
            }
            return std::pair{key, TagValue{valueStr[0]}};
        case 'i': {
            std::int64_t v{0};
            if (!ParseExact(valueStr, v)) {
                return std::nullopt;
            }
            // The spec range for a SAM 'i' tag is [INT32_MIN, UINT32_MAX]; htslib rejects
            // anything outside it rather than truncating into 32 bits.
            if ((v < std::numeric_limits<std::int32_t>::min()) ||
                (v > std::numeric_limits<std::uint32_t>::max())) {
                return std::nullopt;
            }
            return std::pair{key, TagValue{v}};
        }
        case 'f': {
            float v{0.0f};
            if (!ParseExact(valueStr, v)) {
                return std::nullopt;
            }
            return std::pair{key, TagValue{v}};
        }
        case 'Z':
            // SAMv1 §1.5 type 'Z': value must match [ !-~]* (space or printable)
            for (const char ch : valueStr) {
                const auto uch{static_cast<unsigned char>(ch)};
                if (uch != 0x20U && (uch < 0x21U || uch > 0x7EU)) {
                    return std::nullopt;
                }
            }
            return std::pair{key, TagValue{std::string{valueStr}}};
        case 'H':
            if ((std::size(valueStr) % 2) != 0) {
                return std::nullopt;
            }
            if (!std::ranges::all_of(valueStr, IsHexDigit)) {
                return std::nullopt;
            }
            return std::pair{key, TagValue{HexString{std::string{valueStr}}}};
        case 'B': {
            // Format: B:T,v1,v2,...
            if (std::size(valueStr) < 3) {
                return std::nullopt;
            }
            const char elemType{valueStr[0]};
            if (TagArrayElementSize(elemType) == 0) {
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
                std::size_t endPos{std::size(valueStr)};
                if (commaPos != std::string_view::npos) {
                    endPos = commaPos;
                }
                const std::string_view elem{valueStr.substr(pos, endPos - pos)};
                if (elem.empty()) {
                    return std::nullopt;
                }

                if (elemType == 'f') {
                    float v{0.0f};
                    if (!ParseExact(elem, v)) {
                        return std::nullopt;
                    }
                    arr.AppendFloat(v);
                } else {
                    // Integer element types
                    std::int64_t v{0};
                    if (!ParseExact(elem, v)) {
                        return std::nullopt;
                    }

                    if (!AppendIntegerArrayElement(arr, elemType, v)) {
                        return std::nullopt;
                    }
                }

                if (commaPos == std::string_view::npos) {
                    pos = endPos;
                } else {
                    pos = endPos + 1;
                }
            }

            return std::pair{key, TagValue{std::move(arr)}};
        }
        default:
            return std::nullopt;
    }
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
    const std::to_chars_result toCharsResult{
        std::to_chars(std::data(buf), std::data(buf) + std::size(buf), v)};
    out.append(std::data(buf), toCharsResult.ptr);
}

void AppendFloatRaw(std::string& out, float v)
{
    std::array<char, 32> buf{};
    const int n{std::snprintf(std::data(buf), std::size(buf), "%g", static_cast<double>(v))};
    if (n > 0) {
        const std::size_t len{static_cast<std::size_t>(n)};
        out.append(std::data(buf), len);
    }
}

char* AppendUnsignedDecimal(char* dest, std::uint32_t value)
{
    if (value >= 100U) {
        *dest = static_cast<char>('0' + (value / 100U));
        ++dest;
        *dest = static_cast<char>('0' + ((value / 10U) % 10U));
        ++dest;
        *dest = static_cast<char>('0' + (value % 10U));
        ++dest;
        return dest;
    }
    if (value >= 10U) {
        *dest = static_cast<char>('0' + (value / 10U));
        ++dest;
        *dest = static_cast<char>('0' + (value % 10U));
        ++dest;
        return dest;
    }
    *dest = static_cast<char>('0' + value);
    ++dest;
    return dest;
}

char* AppendSignedDecimal(char* dest, std::int32_t value)
{
    if (value < 0) {
        *dest = '-';
        ++dest;
        return AppendUnsignedDecimal(dest, static_cast<std::uint32_t>(-value));
    }
    return AppendUnsignedDecimal(dest, static_cast<std::uint32_t>(value));
}

struct UInt8ArrayOverwriteWriter
{
    const std::byte* data;
    std::uint32_t count;
    std::size_t startPos;

    std::size_t operator()(char* buf, std::size_t bufSize) const;
};

std::size_t UInt8ArrayOverwriteWriter::operator()(char* buf, std::size_t /*bufSize*/) const
{
    char* dest{buf + startPos};
    for (std::uint32_t i{0}; i < count; ++i) {
        *dest = ',';
        ++dest;
        dest = AppendUnsignedDecimal(
            dest, static_cast<std::uint32_t>(static_cast<std::uint8_t>(data[i])));
    }
    return static_cast<std::size_t>(dest - buf);
}

struct Int8ArrayOverwriteWriter
{
    const std::byte* data;
    std::uint32_t count;
    std::size_t startPos;

    std::size_t operator()(char* buf, std::size_t bufSize) const;
};

std::size_t Int8ArrayOverwriteWriter::operator()(char* buf, std::size_t /*bufSize*/) const
{
    char* dest{buf + startPos};
    for (std::uint32_t i{0}; i < count; ++i) {
        const std::int32_t value{static_cast<std::int8_t>(data[i])};
        *dest = ',';
        ++dest;
        dest = AppendSignedDecimal(dest, value);
    }
    return static_cast<std::size_t>(dest - buf);
}

/// \brief Fast-path: serialize B:C array using resize_and_overwrite (no
/// zero-fill).
void SerializeBArrayUInt8(const std::byte* data, std::uint32_t count, std::string& out)
{
    const std::size_t startPos{std::size(out)};
    const std::size_t maxChars{count * 4U};

    out.resize_and_overwrite(startPos + maxChars, UInt8ArrayOverwriteWriter{data, count, startPos});
}

/// \brief Fast-path: serialize B:c array using resize_and_overwrite (no
/// zero-fill).
void SerializeBArrayInt8(const std::byte* data, std::uint32_t count, std::string& out)
{
    const std::size_t startPos{std::size(out)};
    const std::size_t maxChars{count * 5U};

    out.resize_and_overwrite(startPos + maxChars, Int8ArrayOverwriteWriter{data, count, startPos});
}

/// \brief Serialize B:S/s/I/i array elements using to_chars with correct type
/// width.
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
        const std::to_chars_result toCharsResult{
            std::to_chars(std::data(buf), std::data(buf) + std::size(buf), v)};
        out.append(std::data(buf), toCharsResult.ptr);
    }
}

void SerializeBArrayFloat(const std::byte* data, std::uint32_t count, std::string& out)
{
    for (std::uint32_t i{0}; i < count; ++i) {
        out += ',';
        AppendFloatRaw(out, ReadF32LE(data + i * 4));
    }
}

std::size_t AppendNulTerminatedTextTag(std::string& out, char type, std::span<const std::byte> data,
                                       std::size_t offset)
{
    const std::size_t start{offset};
    while ((offset < std::size(data)) && (data[offset] != std::byte{0})) {
        ++offset;
    }
    if (offset == std::size(data)) {
        throw std::runtime_error{"Tags: unterminated BAM auxiliary text payload"};
    }
    const std::span<const std::byte> payload{data.subspan(start, offset - start)};
    if ((type == 'H') &&
        (((std::size(payload) % 2) != 0) || !std::ranges::all_of(payload, [](std::byte value) {
             return IsHexDigit(static_cast<char>(value));
         }))) {
        throw std::runtime_error{"Tags: type 'H' payload is not valid hex"};
    }

    out += type;
    out += ':';
    for (const std::byte value : payload) {
        out += static_cast<char>(value);
    }
    return offset + 1;
}

}  // namespace

void SerializeRawTagsToSam(std::span<const std::byte> data, std::string& out)
{
    std::size_t offset{0};

    while (offset < std::size(data)) {
        RequireBamAuxBytes(data, offset, 3);
        const char c1{static_cast<char>(data[offset])};
        const char c2{static_cast<char>(data[offset + 1])};
        const char type{static_cast<char>(data[offset + 2])};
        offset += 3;

        out += '\t';
        out += c1;
        out += c2;
        out += ':';

        if (type == 'A') {
            RequireBamAuxBytes(data, offset, 1);
            out += "A:";
            out += static_cast<char>(data[offset]);
            offset += 1;
        } else if (type == 'c') {
            RequireBamAuxBytes(data, offset, 1);
            out += "i:";
            AppendIntRaw(out, static_cast<std::int8_t>(data[offset]));
            offset += 1;
        } else if (type == 'C') {
            RequireBamAuxBytes(data, offset, 1);
            out += "i:";
            AppendIntRaw(out, static_cast<std::uint8_t>(data[offset]));
            offset += 1;
        } else if (type == 's') {
            RequireBamAuxBytes(data, offset, 2);
            out += "i:";
            AppendIntRaw(out, ReadI16LE(std::data(data) + offset));
            offset += 2;
        } else if (type == 'S') {
            RequireBamAuxBytes(data, offset, 2);
            out += "i:";
            AppendIntRaw(out, ReadU16LE(std::data(data) + offset));
            offset += 2;
        } else if (type == 'i') {
            RequireBamAuxBytes(data, offset, 4);
            out += "i:";
            AppendIntRaw(out, ReadI32LE(std::data(data) + offset));
            offset += 4;
        } else if (type == 'I') {
            RequireBamAuxBytes(data, offset, 4);
            out += "i:";
            AppendIntRaw(out, ReadU32LE(std::data(data) + offset));
            offset += 4;
        } else if (type == 'f') {
            RequireBamAuxBytes(data, offset, 4);
            out += "f:";
            AppendFloatRaw(out, ReadF32LE(std::data(data) + offset));
            offset += 4;
        } else if (type == 'Z') {
            offset = AppendNulTerminatedTextTag(out, type, data, offset);
        } else if (type == 'H') {
            offset = AppendNulTerminatedTextTag(out, type, data, offset);
        } else if (type == 'B') {
            RequireBamAuxBytes(data, offset, 5);
            const char elemType{static_cast<char>(data[offset])};
            offset += 1;
            const std::uint32_t count{ReadU32LE(std::data(data) + offset)};
            offset += 4;

            out += "B:";
            out += elemType;

            // Fast-path for byte arrays (dominant in PacBio kinetics data)
            if (elemType == 'C') {
                RequireBamAuxBytes(data, offset, count);
                SerializeBArrayUInt8(std::data(data) + offset, count, out);
                offset += count;
            } else if (elemType == 'c') {
                RequireBamAuxBytes(data, offset, count);
                SerializeBArrayInt8(std::data(data) + offset, count, out);
                offset += count;
            } else if (elemType == 's') {
                RequireBamAuxBytes(data, offset, static_cast<std::size_t>(count) * 2);
                SerializeBArrayInt<std::int16_t>(std::data(data) + offset, count, out);
                offset += static_cast<std::size_t>(count) * 2;
            } else if (elemType == 'S') {
                RequireBamAuxBytes(data, offset, static_cast<std::size_t>(count) * 2);
                SerializeBArrayInt<std::uint16_t>(std::data(data) + offset, count, out);
                offset += static_cast<std::size_t>(count) * 2;
            } else if (elemType == 'i') {
                RequireBamAuxBytes(data, offset, static_cast<std::size_t>(count) * 4);
                SerializeBArrayInt<std::int32_t>(std::data(data) + offset, count, out);
                offset += static_cast<std::size_t>(count) * 4;
            } else if (elemType == 'I') {
                RequireBamAuxBytes(data, offset, static_cast<std::size_t>(count) * 4);
                SerializeBArrayInt<std::uint32_t>(std::data(data) + offset, count, out);
                offset += static_cast<std::size_t>(count) * 4;
            } else if (elemType == 'f') {
                RequireBamAuxBytes(data, offset, static_cast<std::size_t>(count) * 4);
                SerializeBArrayFloat(std::data(data) + offset, count, out);
                offset += static_cast<std::size_t>(count) * 4;
            } else {
                throw std::runtime_error{
                    std::format("Tags: unsupported B-array element type '{}'", elemType)};
            }
        } else {
            throw std::runtime_error{std::format("Tags: unsupported BAM tag type '{}'", type)};
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

    void operator()(char v) const
    {
        AppendRawByte(dest, static_cast<std::byte>('A'));
        AppendRawByte(dest, static_cast<std::byte>(v));
    }

    void operator()(std::int64_t v) const
    {
        const char bamType{SmallestIntType(v)};
        AppendRawByte(dest, static_cast<std::byte>(bamType));
        AppendIntegerTagValue(dest, bamType, v);
    }

    void operator()(float v) const
    {
        AppendRawByte(dest, static_cast<std::byte>('f'));
        WriteF32Raw(dest, v);
    }

    void operator()(std::string_view v) const
    {
        AppendRawByte(dest, static_cast<std::byte>('Z'));
        for (const char ch : v) {
            AppendRawByte(dest, static_cast<std::byte>(ch));
        }
        AppendRawByte(dest, std::byte{0});
    }

    void operator()(const HexString& v) const
    {
        AppendRawByte(dest, static_cast<std::byte>('H'));
        for (const char ch : v.value) {
            AppendRawByte(dest, static_cast<std::byte>(ch));
        }
        AppendRawByte(dest, std::byte{0});
    }

    void operator()(const TagArray& v) const
    {
        AppendRawByte(dest, static_cast<std::byte>('B'));
        AppendRawByte(dest, static_cast<std::byte>(v.ElementType()));
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
    for (const auto& [key, value] : tags.Entries()) {
        AppendRawByte(dest, static_cast<std::byte>(key.First()));
        AppendRawByte(dest, static_cast<std::byte>(key.Second()));
        std::visit(BamAppendVisitor{dest}, value);
    }
}

}  // namespace Samoa
}  // namespace PacBio
