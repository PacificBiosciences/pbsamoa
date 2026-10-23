<div align="center">
  <img src="docs/imgs/samoa_logo.png" alt="samoa logo" width="250px"/>
  <h1><i>pbsamoa</i> — SAM Open Alternative</h1>
  <p>A C++23 library for reading and writing SAM/BAM/CRAM files.</p>
  <p>Continue to <a href="docs/README.md">documentation</a> of this early draft.</p>
</div>

## Overview

Built on a decade of bioinformatics experience, _pbsamoa_ is inspired by the work of the _SAM/BAM
Format Specification Working Group_ and _htslib_. It focuses on pipelined I/O for maximum
throughput, deterministic chunking of inputs, and flexible record filtering. Native CRAM v3.x
reader/writer support includes CRAI interop and configurable block/data-series compression methods.
The "open" reflects our commitment to developing this library in the open.

## Highlight: faster sorted BAM merges

`pbsamoa merge` is much faster than `samtools merge` on sorted inputs:

```sh
pbsamoa merge --threads 32 --memory 5G out.bam inputs/*.bam
```

On a 16-input sorted benchmark with 8M HiFi reads:

|       Tool       | Version  | Threads | Run time | Speedup |
| ---------------- | -------- | ------: | -------: | ------: |
| `samtools merge` | 1.23     |      20 |    8m28s |    1.0x |
| `pbsamoa merge`  | e26bab59 |      20 |    2m16s |   ~3.7x |
| `pbsamoa merge`  | e26bab59 |      32 |    1m03s |   ~8.0x |

See [CLI Tools](docs/cli.md#merge--merge-or-concatenate-bam-files) for merge
modes and tuning knobs.

## Quick Start

Everything in this library is in the `PacBio::Samoa` namespace, for clarity prepend all examples
with:
```cpp
using namespace PacBio::Samoa;
```

### `RawRecordView`, `RawRecord`, and `BamRecord`

_pbsamoa_ provides three record types with different trade-offs:

- **`RawRecordView`** borrows raw BAM bytes and decodes fields on demand. Use
  `RawRecordBatch::View()` while the batch remains alive, because `View()` does not copy record bytes.
- **`RawRecord`** owns the raw BAM bytes and decodes fields on demand. Ideal for high-throughput
  single-record APIs or whenever the record must outlive its source buffer.
- **`BamRecord`** stores pre-decoded C++ types (`std::string`, `TagMap`, etc.). Use it when you need
  structured access to many fields, mutation, clipping, or writing records back out.

### BamRawReader – decode-on-demand iteration

`BamRawReader` yields owning `RawRecord` references with no up-front field decoding.
Its batch API provides non-owning `RawRecordView`s without per-record byte copies.

```cpp
#include <pbsamoa/io/BamRawReader.hpp>

#include <print>

BamRawReader reader{"input.bam"};

int mapped{0};
for (const auto& rec : reader.Records()) { mapped += rec.IsMapped(); }

std::println("Mapped reads: {}", mapped);
```

Enable parallel BGZF decompression for large files:

```cpp
#include <pbsamoa/io/BamRawReader.hpp>

BamRawReader reader{"input.bam", BamRawReaderConfig{.BgzfWorkers = 4}};
```

### BamRecordReader – pre-decoded records

`BamRecordReader` wraps `BamRawReader` and decodes every record into a `BamRecord` via a background
thread pool. Fields like `Sequence()` and `Tags()` are ready to use without per-access decoding
cost. The decompressed bytes are omitted.

```cpp
#include <pbsamoa/io/BamRecordReader.hpp>

#include <print>

BamRecordReader reader{"input.bam"};

for (const auto& rec : reader.Records()) {
    if (rec.IsMapped()) {
        std::println("{}\t{}\t{}\t{}", rec.Name(), rec.RefId(), rec.Pos(), rec.MapQ());
    }
}
```

Drop large PacBio pulse tags you don't need:

```cpp
#include <pbsamoa/io/BamRecordReader.hpp>

BamRecordReader reader{"input.bam",
                       BamRecordReaderConfig{.TagFilter = DropTags{TagKey{'i', 'p'}, TagKey{'p', 'w'}}}};
```

### Chunking – parallel-by-ZMW sharding

Chunking splits a BAM file into N disjoint slices by unique `(rgId, zmw)` pairs, so multiple
processes can work on the same file without overlap. It requires a `.zmi` (new) or `.pbi` (legacy)
index alongside the BAM.

```cpp
#include <pbsamoa/io/BamRawReader.hpp>

BamRawReaderConfig cfg{
  .ChunkNum = 2,      // 1-based: this process handles chunk 2
  .TotalChunks = 5,   // divide ZMWs into 5 roughly-equal groups
};
BamRawReader reader{"input.bam", cfg};

for (const auto& rec : reader.Records()) {
    // only records belonging to chunk 2's ZMWs
}
```

**How it works.** The reader opens the ZMW index and counts the unique ZMWs in file order. It divides
them into `TotalChunks` slices with integer boundaries: boundary `k` is `floor(N * k / TotalChunks)`
for `N` ZMWs. Then it seeks the BGZF stream to the first record of the slice and sets a record limit
to the exact count. Reading stops automatically when the chunk is exhausted – no per-record index
lookups after setup.

```
ZMWs:    [10, 20, 30, 40, 50]    5 unique ZMWs, TotalChunks=3
chunk 1: [10]       ← floor(5 * 0 / 3) .. floor(5 * 1 / 3) = indices 0..1
chunk 2: [20, 30]   ← indices 1..3
chunk 3: [40, 50]   ← indices 3..5
```

This is the default `ChunkMode::CONTIGUOUS`. With `ChunkingMode = ChunkMode::SCATTER`, each chunk
gets shuffled tiles of `ChunkTileZmws` ZMWs from all of the file. Scatter mode ignores
`BgzfWorkers` and reads on the caller thread.

Chunking and the config field `Whitelist` are mutually exclusive – the constructor rejects combining
both. The config field does not filter records; use `BamRawReader::Whitelist()` to select ZMWs.
Chunking and `Whitelist` work on one BAM file only. A reader on more than one file rejects them.

### Read Pipeline – layered stages for maximum throughput

_pbsamoa_ composes three pipeline layers. Each layer runs its own threads and hands data to the next
via lock-free SPSC queues, so I/O, decompression, parsing, and decoding all overlap.

```
 Disk
  │  IO thread – reads 32 BGZF blocks per batch
  ▼
 ThreadPool (BgzfWorkers)
  │  libdeflate decompression, one batch per task
  ▼
 Consumer thread – parses raw BAM records from decompressed bytes
  │
  ▼  SPSCQueue<RawRecordBatch>  [8]   one batch = 32 blocks, about 2 MiB
  │
 BamRawReader::ReadBatch(ByteLimit)
  │  gives one RawRecordBatch
  ▼
 Producer thread – reads batches up to BatchBudget,
  │  dispatches decode tasks to ThreadPool (DecodeWorkers)
  │  each record: raw bytes → BamRecord (field decode + tag filter)
  ▼  SPSCQueue<BamRecord>  [OutputCapacity, default 1024]
  │
 Caller thread – BamRecordReader::ReadRecord()
```

**Layer 1 – BGZF (BgzfReader).** Three concurrent stages when `BgzfWorkers > 0`: a dedicated IO
thread reads compressed blocks from disk, a thread pool decompresses them with libdeflate, and a
consumer thread parses the decompressed bytes into one `RawRecordBatch` for each group of 32 blocks.
With `BgzfWorkers = 0` everything runs synchronously on the caller thread.

**Layer 2 – Batching (BamRawReader).** `ReadBatch(ByteLimit)` gives a `RawRecordBatch` – one
record-byte buffer plus N extents. The default limit is 256 MiB. With `BgzfWorkers = 0`, the batch
copies records from Layer 1 until it reaches the limit. With `BgzfWorkers > 0`, a batch holds the
records of at most one group of 32 blocks (about 2 MiB). A larger limit does not join groups.
`RawRecordBatch::View()` decodes directly from that buffer and does not allocate a byte copy. This
amortises allocation overhead and gives Layer 3 a batch-sized unit of work.

**Layer 3 – Decoding (BamRecordReader).** A background producer thread calls `ReadBatch` until it
has `BatchBudget` bytes (default 4 MiB). Then it fans out record-level decoding across
`DecodeWorkers` (default 4) threads. Each record is independently converted from raw BAM bytes to a
`BamRecord` with structured fields, applying `TagFilter` (drop/keep) during `ToOwned()`. Decoded
records are pushed into the output SPSC queue for the caller.

Every queue boundary has a named stall counter, so you can identify the bottleneck stage without
external instrumentation. `BamRawReader::GetMetrics()` returns `BgzfMetrics` with `IoStalls`,
`ConsumerStalls`, and `ReaderStalls`. `BamRecordReader::GetMetrics()` returns `ReaderMetrics`. It
adds `DecodeMetrics` with `ProducerStalls` and `ConsumerStalls`. With `BgzfWorkers = 0`, only
`RecordsConsumed` of `BgzfMetrics` has a value.

```cpp
#include <pbsamoa/io/BamRecordReader.hpp>

using namespace PacBio::Samoa::Literals;

BamRecordReaderConfig cfg;
cfg.RawReaderConfig.BgzfWorkers = 8;  // layer 1: decompression threads
cfg.BatchBudget = 4_MiB;              // layer 3: bytes per decode round
cfg.DecodeWorkers = 4;                // layer 3: decode threads
cfg.TagFilter = DropTags{TagKey{'i', 'p'}, TagKey{'p', 'w'}};
BamRecordReader reader{"input.bam", cfg};
```
