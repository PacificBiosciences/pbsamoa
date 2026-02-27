<div align="center">
  <img src="docs/imgs/samoa_logo.png" alt="samoa logo" width="250px"/>
  <h1>pbsamoa : SAM Open Alternative</h1>
  <p>A C++23 library for reading and writing SAM/BAM files.</p>
  </p>Continue to <a href="docs/README.md">documentation</a> of this early draft.</p>
</div>

## Quick Start

### RawRecord vs BamRecord

pbsamoa provides two record types with different trade-offs:

- **`RawRecord`** owns the raw BAM bytes and decodes fields on demand. Ideal for
  high-throughput scanning where you only inspect a few fields per record
  (filtering, counting, indexing).
- **`BamRecord`** stores pre-decoded C++ types (`std::string`, `TagMap`, etc.).
  Use it when you need structured access to many fields, mutation, clipping, or
  writing records back out.

### BamRawReader – zero-copy iteration

`BamRawReader` yields `RawRecord` references directly from decompressed BGZF
blocks with no up-front decoding.

```cpp
#include <pbsamoa/io/BamRawReader.hpp>

#include <print>

int main()
{
    PacBio::Samoa::BamRawReader reader{"input.bam"};

    int mapped{0};
    for (const auto& rec : reader.Records()) {
        if (rec.IsMapped()) {
            ++mapped;
        }
    }

    std::println("Mapped reads: {}", mapped);
}
```

Enable parallel BGZF decompression for large files:

```cpp
PacBio::Samoa::BamRawReaderConfig cfg;
cfg.BgzfWorkers = 4;
PacBio::Samoa::BamRawReader reader{"input.bam", cfg};
```

### BamRecordReader – pre-decoded records

`BamRecordReader` wraps `BamRawReader` and decodes every record into a
`BamRecord` via a background thread pool. Fields like `Sequence()` and `Tags()`
are ready to use without per-access decoding cost.

```cpp
#include <pbsamoa/io/BamRecordReader.hpp>

#include <print>

int main()
{
    PacBio::Samoa::BamRecordReader reader{"input.bam"};

    for (const auto& rec : reader.Records()) {
        if (rec.IsMapped()) {
            std::println("{}\t{}\t{}\t{}",
                         rec.Name(), rec.RefId(), rec.Pos(), rec.MapQ());
        }
    }
}
```

Drop large PacBio pulse tags you don't need:

```cpp
PacBio::Samoa::BamRecordReaderConfig cfg;
cfg.TagFilter = PacBio::Samoa::DropTags{{"ip", "pw"}};
PacBio::Samoa::BamRecordReader reader{"input.bam", cfg};
```

### Chunking – parallel-by-ZMW sharding

Chunking splits a BAM file into N disjoint slices by unique `(rgId, zmw)` pairs,
so multiple processes can work on the same file without overlap. It requires a
`.zmi` or `.pbi` index alongside the BAM.

```cpp
PacBio::Samoa::BamRawReaderConfig cfg;
cfg.ChunkNum = 2;      // 1-based: this process handles chunk 2
cfg.TotalChunks = 5;   // divide ZMWs into 5 roughly-equal groups
PacBio::Samoa::BamRawReader reader{"input.bam", cfg};

for (const auto& rec : reader.Records()) {
    // only records belonging to chunk 2's ZMWs
}
```

**How it works.** The reader opens the ZMW index, gets the file-order list of unique
ZMWs, divides them into `TotalChunks` slices with floating-point boundaries
(`std::lround`), then seeks the BGZF stream to the first record of the slice
and sets a record limit to the exact count. Reading stops automatically when the
chunk is exhausted – no per-record index lookups after setup.

```
ZMWs:    [10, 20, 30, 40, 50]    5 unique ZMWs, TotalChunks=3
chunk 1: [10, 20]   ← lround(5/3 * 0) .. lround(5/3 * 1) = indices 0..2
chunk 2: [30, 40]   ← indices 2..3
chunk 3: [50]       ← indices 3..5 (clamped)
```

Chunking and `Whitelist` are mutually exclusive – the constructor rejects
combining both.

### Read Pipeline – layered stages for maximum throughput

pbsamoa composes three pipeline layers. Each layer runs its own threads and
hands data to the next via lock-free SPSC queues, so I/O, decompression,
parsing, and decoding all overlap.

```
 Disk
  │  IO thread – reads 32 BGZF blocks per batch
  ▼
 ThreadPool (BgzfWorkers)
  │  libdeflate decompression, one batch per task
  ▼
 Consumer thread – parses raw BAM records from decompressed bytes
  │
  ▼  SPSCQueue<RawRecord>  [4096]
  │
 BamRawReader::ReadBatch(ByteLimit)
  │  accumulates records into a single contiguous RawRecordBatch
  ▼
 Producer thread – dispatches decode tasks to ThreadPool (DecodeWorkers)
  │  each record: RawRecord → BamRecord (field decode + tag filter)
  ▼  SPSCQueue<BamRecord>  [4096]
  │
 Caller thread – BamRecordReader::ReadRecord()
```

**Layer 1 – BGZF (BgzfReader).** Three concurrent stages when `BgzfWorkers > 0`:
a dedicated IO thread reads compressed blocks from disk, a thread pool
decompresses them with libdeflate, and a consumer thread parses the
decompressed bytes into `RawRecord`s. With `BgzfWorkers = 0` everything runs
synchronously on the caller thread.

**Layer 2 – Batching (BamRawReader).** `ReadBatch(ByteLimit)` drains records
from Layer 1 into a contiguous `RawRecordBatch` – one allocation, N extents –
up to a memory budget (default 256 MiB). This amortises allocation overhead and
gives Layer 3 a batch-sized unit of work.

**Layer 3 – Decoding (BamRecordReader).** A background producer thread calls
`ReadBatch`, then fans out record-level decoding across `DecodeWorkers` (default
4) threads. Each record is independently converted from raw BAM bytes to a
`BamRecord` with structured fields, applying `TagFilter` (drop/keep) during
`ToOwned()`. Decoded records are pushed into the output SPSC queue for the
caller.

Every queue boundary has a named stall counter (`ioStalls`, `consumerStalls`,
`producerStalls`, etc.) accessible from `GetMetrics()`, so you can identify the
bottleneck stage without external instrumentation.

```cpp
PacBio::Samoa::BamRecordReaderConfig cfg;
cfg.BgzfWorkers = 8;      // layer 1: decompression threads
cfg.DecodeWorkers = 4;     // layer 3: decode threads
cfg.BatchBudget = 4_MiB;   // layer 2: batch size
cfg.TagFilter = PacBio::Samoa::DropTags{{"ip", "pw"}};
PacBio::Samoa::BamRecordReader reader{"input.bam", cfg};
```
