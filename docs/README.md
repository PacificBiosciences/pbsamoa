# pbsamoa : SAM Open Alternative

A C++23 library for reading and writing SAM/BAM files, with BAI index
support. Built for throughput: parallel BGZF decompression, zero-copy record
views, and batch processing with memory budgets.

External dependencies: [libdeflate](https://github.com/ebiggers/libdeflate).

## Building

```sh
meson setup build
meson compile -C build
meson test -C build
```

Dependencies are resolved via Meson subproject wraps — no system packages
required.

## Quick Start

### Read a BAM file (views)

```cpp
#include <pbsamoa/io/BamRawReader.hpp>
#include <print>

PacBio::Samoa::BamRawReader reader{"input.bam"};

// Print every read name (zero-copy views)
for (const auto& view : reader.Records()) {
    std::println("{}", view.Name());
}
```

### Read a BAM file (owned records)

```cpp
#include <pbsamoa/io/BamRecordReader.hpp>
#include <print>

PacBio::Samoa::BamRecordReader reader{"input.bam"};

// Pre-decoded owned BamRecord objects via background pipeline
for (const auto& record : reader.Records()) {
    std::println("{}", record.Name());
}
```

### Read in batches (parallel decompression)

```cpp
#include <pbsamoa/io/BamRawReader.hpp>

using namespace PacBio::Samoa::Literals;

PacBio::Samoa::BamRawReader reader{"input.bam",
    PacBio::Samoa::BamRawReaderConfig{.BgzfWorkers = 8}};

while (auto batch = reader.ReadBatch(128_MiB)) {
    for (std::size_t i{0}; i < batch->RecordCount(); ++i) {
        const PacBio::Samoa::RawRecord view{batch->RecordData(i)};
        // process view
    }
}
```

### Write a BAM file

```cpp
#include <pbsamoa/io/BamRawReader.hpp>
#include <pbsamoa/io/BamWriter.hpp>

PacBio::Samoa::BamRawReader reader{"input.bam"};
PacBio::Samoa::BamWriter writer{"output.bam", reader.Header()};

for (const auto& view : reader.Records()) {
    writer.Write(view);  // zero-copy passthrough
}
```

### Region query with BAI index

```cpp
#include <pbsamoa/io/BamRawReader.hpp>
#include <pbsamoa/index/BaiIndex.hpp>

PacBio::Samoa::BamRawReader reader{"sorted.bam"};
auto index = PacBio::Samoa::BaiIndex::FromFile("sorted.bam.bai");

std::int32_t refId{reader.Header().ReferenceId("chr1")};

for (const auto& view : reader.Query(index, refId, 1000, 2000)) {
    // records overlapping chr1:1000-2000
}
```

### Mutate records

```cpp
#include <pbsamoa/io/BamRawReader.hpp>
#include <pbsamoa/io/BamWriter.hpp>

PacBio::Samoa::BamRawReader reader{"input.bam"};
PacBio::Samoa::BamWriter writer{"output.bam", reader.Header()};

for (const auto& view : reader.Records()) {
    auto record = view.ToOwned();
    record.MapQ(60);
    writer.Write(record);
}
```

### Filter tags during conversion

Large tags (PacBio kinetics `ip`/`pw`) can dominate memory. Drop them when
converting views to owned records:

```cpp
auto record = view.ToOwned(DropTags{TagKey{'i', 'p'}, TagKey{'p', 'w'}});
```

Or keep only the tags you need:

```cpp
auto record = view.ToOwned(KeepTags{TagKey{'R', 'G'}, TagKey{'N', 'M'}});
```

## Documentation

|            Document             |                  Contents                  |
| ------------------------------- | ------------------------------------------ |
| [Architecture](architecture.md) | Layer diagram, design decisions, data flow |
| [API Reference](api.md)         | All public types and functions             |
| [CLI Tools](cli.md)             | The `pbsamoa` command-line tool            |
| [C++23 Guide](cpp23-best-practices.md) | Project-level C++23 coding rules and applicability |
