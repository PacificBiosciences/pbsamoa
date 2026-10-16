# pbsamoa : SAM Open Alternative

A C++23 library for reading and writing SAM/BAM/CRAM files, with BAI/CRAI
index support. Built for throughput: parallel BGZF decompression, zero-copy
record views, and batch processing with memory budgets.

Required dependencies:

- [libdeflate](https://github.com/ebiggers/libdeflate) (Meson wrap, or a system package)
- pbcopper (Meson wrap). pbcopper needs Boost from the system and gets simde as a wrap.
- threads

Optional dependencies:

- bzip2 and liblzma (system packages). They enable the BZIP2 and LZMA block methods of CRAM.
- [ISA-L](https://github.com/intel/isa-l) (Meson wrap, x86-64 only, needs NASM). It makes BGZF
  compression faster. Use the option `-Disal=enabled|disabled|auto` (default `auto`).
- pbtestutils and googletest (Meson wraps) for the tests. Use `-Dtests=false` to skip them.

The source tree includes [htscodecs](https://github.com/samtools/htscodecs) and the rigtorp
`SPSCQueue` in `third-party/`.

## Building

```sh
meson setup build
meson compile -C build
meson test -C build
```

Meson gets libdeflate, pbcopper, ISA-L, and the test dependencies through subproject wraps.
Install Boost, and optionally bzip2 and liblzma, as system packages.

## Quick Start

### Read a BAM file (raw records)

```cpp
#include <pbsamoa/io/BamRawReader.hpp>
#include <print>

PacBio::Samoa::BamRawReader reader{"input.bam"};

// Print every read name (decode on demand)
for (const auto& record : reader.Records()) {
    std::println("{}", record.Name());
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

#include <cstddef>

using namespace PacBio::Samoa::Literals;

PacBio::Samoa::BamRawReader reader{"input.bam",
    PacBio::Samoa::BamRawReaderConfig{.BgzfWorkers = 8}};

while (auto batch = reader.ReadBatch(128_MiB)) {
    for (std::size_t i{0}; i < batch->RecordCount(); ++i) {
        const PacBio::Samoa::RawRecordView view{batch->View(i)};
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

for (const auto& record : reader.Records()) {
    writer.Write(record);  // writes the raw bytes, no encode step
}
```

`Records()` makes an owning `RawRecord` for each record. To skip this step, write whole
batches:

```cpp
#include <pbsamoa/io/BamRawReader.hpp>
#include <pbsamoa/io/BamWriter.hpp>

PacBio::Samoa::BamRawReader reader{"input.bam"};
PacBio::Samoa::BamWriter writer{"output.bam", reader.Header()};

while (auto batch = reader.ReadBatch()) {
    writer.WriteBatch(*batch);
}
```

### Read a CRAM file

```cpp
#include <pbsamoa/io/CramReader.hpp>
#include <print>

PacBio::Samoa::CramReader reader{"input.cram"};
for (const auto& record : reader.Records()) {
    std::println("{}", record.Name());
}

// Raw BAM-layout access from CRAM decode
for (const auto& raw : reader.RawRecords()) {
    std::println("{}", raw.Name());
}
```

### Write a CRAM file

```cpp
#include <pbsamoa/io/BamRecordReader.hpp>
#include <pbsamoa/io/CramWriter.hpp>

PacBio::Samoa::BamRecordReader reader{"input.bam"};
PacBio::Samoa::CramWriter writer{
    "output.cram",
    reader.Header(),
    PacBio::Samoa::CramWriterConfig{.WriteCrai = true},
};

for (const auto& record : reader.Records()) {
    writer.Write(record);
}
writer.Close();
```

### Region query with BAI index

```cpp
#include <pbsamoa/index/BaiIndex.hpp>
#include <pbsamoa/io/BamRawReader.hpp>

#include <cstdint>

PacBio::Samoa::BamRawReader reader{"sorted.bam"};
auto index = PacBio::Samoa::BaiIndex::FromFile("sorted.bam.bai");

std::int32_t refId{reader.Header().ReferenceId("chr1")};

for (const auto& record : reader.Query(index, refId, 1000, 2000)) {
    // records that overlap chr1 [1000, 2000), 0-based, end exclusive
}
```

### Mutate records

```cpp
#include <pbsamoa/io/BamRawReader.hpp>
#include <pbsamoa/io/BamWriter.hpp>

PacBio::Samoa::BamRawReader reader{"input.bam"};
PacBio::Samoa::BamWriter writer{"output.bam", reader.Header()};

for (const auto& raw : reader.Records()) {
    auto record = raw.ToOwned();
    record.MapQ(60);
    writer.Write(record);
}
```

### Filter tags during conversion

Large tags (PacBio kinetics `ip`/`pw`) can dominate memory. Drop them when
converting views to owned records:

```cpp
#include <pbsamoa/core/RawRecord.hpp>
#include <pbsamoa/core/Tags.hpp>

using namespace PacBio::Samoa;

auto record = view.ToOwned(DropTags{TagKey{'i', 'p'}, TagKey{'p', 'w'}});
```

Or keep only the tags you need:

```cpp
#include <pbsamoa/core/RawRecord.hpp>
#include <pbsamoa/core/Tags.hpp>

using namespace PacBio::Samoa;

auto record = view.ToOwned(KeepTags{TagKey{'R', 'G'}, TagKey{'N', 'M'}});
```

## Documentation

|            Document             |                  Contents                  |
| ------------------------------- | ------------------------------------------ |
| [Architecture](architecture.md) | Layer diagram, design decisions, data flow |
| [API Reference](api.md)         | All public types and functions             |
| [CLI Tools](cli.md)             | The `pbsamoa` command-line tool            |
