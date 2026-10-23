# CLI Tools

pbsamoa ships a multi-tool binary with subcommands for BAM/SAM/CRAM workflows,
indexing, and performance benchmarking.

## pbsamoa

```sh
pbsamoa <tool> [options] <arguments>
pbsamoa --help              # list the tools
pbsamoa <tool> --help       # options of one tool
pbsamoa --version
```

The tools are `dump`, `chunk`, `sort`, `merge`, `convert`, `subsample`,
`depth`, `undo`, `bai-build`, `bai-query`, `zmi-build`, `zmi-index`,
`zmi-query`, and `bench`.

### Common options and exit status

All tools have these options:

| Flag | Description |
| ---- | ----------- |
| `-h`, `--help` | Show the help and exit |
| `--version` | Show the version and exit |
| `--log-level LEVEL` | `TRACE`, `DEBUG`, `INFO`, `WARN`, or `FATAL` (default: `WARN`) |
| `--log-file FILE` | Write the log to `FILE`, not to stderr |

Tool options have only long names. `-h` is the only short flag.

A tool exits with status `0` on success. On an error, it writes the message
to stderr and exits with status `1`. An unknown tool name, an unknown option,
or too many arguments also give status `1`. A tool run without its arguments
(`pbsamoa sort`), and `pbsamoa` alone, print the help and exit with status `0`.

Memory sizes (`--memory`, `--batch-bytes`) accept a `K`, `M`, or `G` suffix
in upper or lower case. These suffixes are powers of 1024. The `--yield` of
`subsample` is different: its suffixes are powers of 1000.

### dump — BAM/SAM/CRAM to SAM text

Convert BAM, SAM, or CRAM input to SAM text on stdout.

```sh
pbsamoa dump input.bam
pbsamoa dump input.sam
pbsamoa dump input.cram
pbsamoa dump input.cram > output.sam
pbsamoa dump --header-only input.bam
```

Options:

| Flag | Description |
| ---- | ----------- |
| `--bgzf-threads N` | Decompression worker threads; `-1` = auto, `min(hw, 8)` (default: `-1`). An explicit count is capped at 10 |
| `--format-threads N` | SAM formatting thread pool size for BAM input; `0` = auto, `max(hw, 4)` (default: `0`) |
| `--reference ref.fa` | Reference FASTA for reference-based CRAM decoding. It is memory-mapped and read where needed, at offsets from `ref.fa.fai` (without one, the file is scanned once) |
| `--region REGION` | CRAM-only region query: `ref`, `ref:start`, or `ref:start-end` (1-based, inclusive), or `*` for the unmapped records |
| `--index file.crai` | CRAM-only CRAI path override (requires `--region`) |
| `--no-header` | Do not write the SAM header |
| `--header-only` | Write only the SAM header |

Notes:

- `--region` and `--index` are only supported for CRAM input.
- If `--region` is set and `--index` is omitted, `dump` loads `<input>.crai`.
- A region splits on its rightmost colon, so reference names may contain `:`.
  `ref` is the whole reference and `ref:start` runs to its end.
- `--no-header` and `--header-only` are mutually exclusive.
- The tool selects the input format from the file extension: `.sam`, `.cram`,
  or else BAM.

Environment:

- `PBSAMOA_METRICS=1` — Print per-second pipeline metrics and a final summary to stderr for BAM dump pipeline.

### convert — BAM/SAM to CRAM

Convert BAM or SAM input to CRAM output.

```sh
pbsamoa convert input.bam output.cram
pbsamoa convert input.sam output.cram
pbsamoa convert --reference ref.fa input.bam output.cram   # reference-based
```

Options:

| Flag | Description |
| ---- | ----------- |
| `--records-per-slice N` | Maximum records per CRAM slice (default: `10000`; must be > 0) |
| `--bgzf-threads N` | BAM input decompression workers (default: `-1` = auto, max 8; `0` = serial; explicit counts capped at 16) |
| `--convert-to-bam-record` | Opt in to decoding BAM input into `BamRecord` before CRAM write |
| `--decode-threads N` | BAM decode workers for `--convert-to-bam-record` mode (default: `-1` = auto, max 8; `0` = serial; explicit counts capped at 16). A value of `0` or more also enables that mode |
| `--compression-threads N` | CRAM block compression workers (default: `-1` = auto, max 8; `0` = serial; explicit counts capped at 8) |
| `--block-compression METHOD` | Default CRAM block method (default: `rans4x8`) |
| `--series-compression SERIES=METHOD` | Override method for a CRAM data series (repeatable; `SERIES` is not case-sensitive) |
| `--write-crai` | Write samtools-compatible `<output>.crai` |
| `--reference ref.fa` | Write reference-based CRAM: mapped records store only their differences to the reference, and @SQ lines without `M5` get one (an existing `M5` is kept). Readers then need the reference. A file without sequences (empty, or gzip-compressed) is an error. The FASTA is memory-mapped; with `ref.fa.fai` nothing but the needed bases and the `M5` hashes reads it |

Supported `METHOD` values:

- `raw`, `gzip`, `bzip2`, `lzma`, `rans4x8`, `rans4x16`, `arith`, `fqzcomp`, `tok`
- Numeric method IDs `0` through `8`

Supported `SERIES` values:

- `BF`, `CF`, `RI`, `RL`, `AP`, `RG`, `RN`, `MF`, `NS`, `NP`, `TS`, `NF`, `TL`, `FN`, `FC`, `FP`, `DL`, `BB`, `QQ`, `BS`, `IN`, `RS`, `PD`, `HC`, `SC`, `MQ`, `BA`, `QS`

Notes:

- Output must use `.cram` extension.
- Input must use `.bam` or `.sam` extension.
- Default block compression method is `rans4x8`.
- BAM input defaults to `RawRecord` passthrough into `CramWriter`; use
  `--convert-to-bam-record` to use the decoded `BamRecord` path.
- Advanced `CramWriterConfig` fields (`SlicesPerContainer`,
  `CompressionLevel`, `UseTempFile`) are currently API-only and are not
  exposed as `convert` flags.

### undo — undo barcode clipping

Reverse demultiplexing done by lima or horizon: re-attach the clipped barcode
flanks (sequence, qualities, kinetics, MM/ML, subread pileup) from the `ls` undo
stack, unwinding nested clips, and drop the demultiplexing tags
(`bc bq bx bl bt ql qt ls`, CCS `qs`/`qe`, barcode bits in `cx`). Barcoded read
groups collapse onto their base ID and lima `@PG` records are removed. A record
without an `ls` stack is not restored: only its `RG` collapses, and it keeps its
other demultiplexing tags.

```sh
pbsamoa undo demux.bam original.bam
```

The tool has no options other than the common options. It adds a
`pbsamoa.undo` `@PG` record. On stderr it prints
`pbsamoa undo: N records, M restored`.

### subsample — keep a fraction of the templates, or a target yield

Keep a fraction of the templates in a BAM file, or keep random ZMWs until their
bases reach a target yield. Give `--fraction` or `--yield`, not both.

```sh
pbsamoa subsample --fraction 0.1 input.bam output.bam            # auto seed
pbsamoa subsample --fraction 0.1 --seed 42 input.bam output.bam  # fixed seed
pbsamoa subsample --yield 60G hifi_reads.bam output.bam          # 60 * 10^9 bases
```

Options:

| Flag | Description |
|------|-------------|
| `--fraction F` | Fraction of templates to keep, in (0, 1] |
| `--yield Y` | Target yield in bases; suffixes K, M, G, T are powers of 1000 |
| `--seed S` | Unsigned 32-bit integer, or `auto` (default): the FNV-1a hash of the input header text, as in samtools >= 1.24 |
| `--threads N` | BGZF threads for reading and for writing; -1 = auto, `min(hw, 8)` (default). An explicit count is capped at 64 |
| `--compression N` | Output BGZF compression level in [1, 12] (default 6) |
| `--strip-kinetics` | Remove the HiFi kinetics tags `fi`, `fp`, `ri`, `rp`, `ip`, and `pw` from the kept records |
| `--status-json PATH` | Write a live status stream to PATH (see below) |

Both modes write a `@PG` record. The tool prints the kept records and bases on
stderr. On an error, the tool removes the partial output BAM.

**`--status-json PATH`** writes one JSON object per line: `run_info` first, a
`start` and an `end` line for each phase, a `progress` line each second while
the count changes, and `final` last (`FINISHED`, or `FAILED` with the cause).
`run_info` contains the command line, the host name, and the thread count.
The phases:

| Mode | Phases | Unit and total of `progress` |
|------|--------|------------------------------|
| `--yield` with a `.pbi` | `select`, `write` | `select`: none; `write`: reads, total from the `.pbi` |
| `--yield` without a `.pbi` | `count`, `write` | `count`: reads, total not known; `write`: reads, total from `count` |
| `--fraction` | `write` | reads, total from the `.pbi` if there is one |

The `end` lines have the counts as `metrics`. The `end` line of `select` or
`count` has `records`, `bases`, and `target_bases`. The `end` line of `write`
has `records`, `bases`, `kept_records`, and `kept_bases`. `final` has the
`write` counts, `target_bases` or `fraction`, `seed`, and `strip_kinetics`
when it is set. `final` lists the output BAM in `outputs`. A write error disables the stream
with one warning; it never fails the run.

**`--fraction`** replaces `samtools view --subsample FRAC --subsample-seed SEED`
and keeps the same records. The tool hashes the read name, so all records of a
template (mates, subreads of one ZMW name) are kept or dropped together. The
output header gets `@CO Sub-sampled fraction=F seed=S`, as samtools writes it.
Like samtools, the tool scrambles a non-zero seed with the C library `rand()`.
glibc and macOS give different values, so for a non-zero seed the kept records
are the same as samtools on the same platform only. Seed 0 gives the same
records on every platform. Unlike samtools, `--fraction 0` is an error, not
"keep everything".

**`--yield`** is for unaligned HiFi reads. The yield is the sum of the read
lengths. The tool hashes (read group, ZMW) with the seed into 2^24 buckets and
keeps the buckets in hash order until their bases reach the target. Thus the
kept yield is at least the target, and it overshoots by less than one bucket,
usually one read. Both strand records of a by-strand ZMW stay together. A
larger target with the same seed keeps a superset of the reads.

- With `BAM.pbi`, the tool reads the lengths from the index (`qEnd - qStart`)
  and reads the BAM once. It checks each record against its index row and stops
  with an error if the `.pbi` is stale. Then it removes the partial output.
- Without a `.pbi`, the tool reads the BAM twice: once for the lengths and once
  to write. The kept records are the same as with the `.pbi`.
- Read names must be PacBio names `movie/zmw/...`.
- The selection does not use `rand()`, so it is the same on every platform. It
  does not match a samtools run, because samtools has no yield mode.
- If the input has less than the target, the tool keeps every record and prints
  a warning.
- The bucket table uses a fixed 128 MiB of memory (2^24 bucket counts of 8
  bytes), for any input size.

The output header gets `@CO Sub-sampled yield=Y seed=S`. `Y` is the number
of bases, for example `yield=60000000000` for `--yield 60G`.

### depth — mean read depth in bins

Write the mean read depth of each bin of each contig. The tool counts the same
positions as `samtools depth -J`: `M`, `=`, `X` and `D` count, `N` does not.
Unmapped, secondary, QC-fail and duplicate records do not count; supplementary
records and MAPQ 0 do.

```sh
pbsamoa depth aligned.bam depth.tsv                  # 1 kb bins
pbsamoa depth --bin-size 10000 aligned.bam depth.tsv
```

Options:

| Flag | Description |
|------|-------------|
| `--bin-size N` | Bin size in bp (default 1000; must be at least 1) |
| `--threads N` | BGZF threads and depth-sum threads; -1 = auto, `min(hw, 8)` (default). An explicit count is capped at 64 |
| `--status-json PATH` | Write a live status stream to PATH, as `subsample` |

The output is a TSV file with a `#contig  start  end  mean_depth` header line.
It has one line for each bin of each contig, in header order, also for bins and
contigs without reads. `start` is 0-based and `end` is exclusive. The last bin
of a contig is shorter. The mean prints as the shortest text that reads back as
the same double, so `mean_depth * (end - start)` gives the exact integer sum of
the depth in the bin.

The tool reads the BAM once and does not need an index or a sorted input. Each
thread keeps one 8-byte sum for each bin: for GRCh38 at 1 kb, that is 25 MB for
each thread. The status stream has one phase, `depth`, in reads. With
`BAM.bai`, the progress has a total: the record count of the index (the mapped
and unmapped counts of all references and the unplaced reads). A `.bai` older
than the BAM gives no total. The metrics are `records`, `counted_records`,
`depth_bases` (the sum of the depth over all positions), `bins`, `bin_size`,
and `index_records` when the `.bai` gives the total. `final` lists the TSV
file in `outputs`. On stderr the tool prints
`pbsamoa depth: N of M records counted, D depth bases in B bins of S bp`.

### sort — sort a BAM file

Sort a BAM file by coordinate, query name, or a 2-char tag, using a RAM-bounded
external merge sort (spill to temp run files, then a bounded-fan-in k-way merge —
multi-pass when the run count exceeds the open-file limit).

```sh
pbsamoa sort input.bam sorted.bam
pbsamoa sort --order queryname input.bam byname.bam
```

Options:

| Flag | Description |
| ---- | ----------- |
| `--order ORDER` | `coordinate` \| `queryname` \| `tag` (default: `coordinate`) |
| `--tag XX` | 2-char tag to sort by (required iff `--order tag`) |
| `--memory SIZE` | Peak run-buffer budget, `K`/`M`/`G` suffix (default: `768M`) |
| `--temp-dir DIR` | Directory for temporary run files (default: output directory) |
| `--threads N` | Worker threads, `0`/auto = `min(hw, 8)` (default: `0`) |
| `--compression L` | Output BGZF level `[1,12]` (default: `6`) |

The tool adds a `pbsamoa.sort` `@PG` record. On stderr it prints
`pbsamoa sort: N records, M spilled run(s)`.

### merge — merge or concatenate BAM files

Combine several BAM files (which must share an identical `@SQ` reference list)
into one. The merged header unions their `@RG` / `@PG` / `@CO` records and an
appended `pbsamoa.merge` `@PG`. With no `--order` / `--concat`, the mode is
auto-detected from the inputs' `@HD SO`:

- all `coordinate` or all `queryname` → streaming k-way sorted merge of the
  inputs (each must already be sorted by that order); output is deterministic
  and independent of the thread count, decode timing, and `--memory` budget.
- all unsorted/unknown → byte-level concatenation.
- a mix → error (pass `--order` or `--concat`).

Performance note: `pbsamoa merge` is much faster than `samtools merge` on
sorted inputs. On a 16-input benchmark with 8,535,828 records, `samtools 1.23
merge -@20` took 8m28s, while `pbsamoa merge --threads 20 --memory 5G` took
2m16s (~3.7x faster). With `--threads 32 --memory 5G`, `pbsamoa merge`
finished in 1m03s (~8.0x faster than that samtools run).

The sorted merge reads every input in parallel and ahead of the merge: each
input is decompressed on a shared worker pool and framed into a per-input bounded
queue, while the heap consumes already-decoded records and the writer compresses
the output in parallel. `--memory` caps the total in-flight read-ahead across all
inputs (a single budget, independent of input count); each input always keeps its
head record available, so a budget below one record cannot stall the merge.

When a coordinate merge's inputs occupy strictly disjoint coordinate ranges — each
file's coordinates entirely precede the next's — the merge skips the heap entirely
and emits via verbatim BGZF block passthrough (the `--concat` machinery), reordered
by minimum coordinate, after `ProbeDisjointChain` verifies non-overlap and internal
sort order. This eliminates recompression — the dominant cost — so the common shard
/ per-chromosome / `chunk`→`merge` roundtrip runs several times faster than the heap
merge (and than `samtools merge`, which always recompresses). Overlapping, internally
unsorted, or query-name/tag merges use the heap path. The CLI summary notes
`(passthrough)` when this path is taken.

```sh
pbsamoa merge out.bam a.bam b.bam              # auto-detect
pbsamoa merge --order coordinate out.bam a.bam b.bam
pbsamoa merge --concat out.bam a.bam b.bam     # force concatenation
```

`--concat` concatenates record streams via BGZF block passthrough — compressed
blocks are copied input→output without decompress/recompress (only the single
block straddling each input's header/record boundary is rebuilt), which is the
fastest path. The output is marked `SO:unsorted` and its record count is not
reported (counting would require decompressing every block).

Options:

| Flag | Description |
| ---- | ----------- |
| `--order ORDER` | `coordinate` \| `queryname` \| `tag` (default: auto-detect) |
| `--concat` | Byte-concatenate inputs (BGZF passthrough), output `SO:unsorted` |
| `--tag XX` | 2-char tag inputs are sorted by (required iff `--order tag`) |
| `--threads N` | Default size of both CPU pools (input decode + output compress), `0`/auto = `min(hw, 8)`. Does **not** bound the per-input producer threads (one I/O-bound thread per input). |
| `--decode-threads N` | Input-decompression pool size; `0` = inherit `--threads`. Sorted merge only. |
| `--compress-threads N` | Output-compression pool size; `0` = inherit `--threads`. Sorted merge only. |
| `--memory SIZE` | Total input read-ahead budget, `K`/`M`/`G` suffix (default: `768M`); sorted merge only. **Not a total memory cap**: peak RSS runs ~1.4x the budget once producers pin it, plus ~`NumInputs × --batch-bytes` (per-input current batches held outside the budget), the writer queue, and the `--bai` builder. Under a hard limit (cgroup/SLURM), set the limit to at least 2x `--memory`. |
| `--batch-bytes SIZE` | Per-input handoff batch target, `K`/`M`/`G` suffix (default: `256K`); sorted merge only. |
| `--writer-queue N` | Output writer queue depth in compressed blocks (default: `256`, must be `≥ 1`); sorted merge only. |
| `--compression L` | Output BGZF level `[1,12]` (default: `6`) |
| `--bai` | Write `<out>.bai` during the merge, built on the fly (no extra pass). Coordinate merges that take the heap path only — rejected (non-zero exit) for `--concat`, `queryname`/`tag` order, or the disjoint-chain passthrough, since none decodes records into an indexable coordinate-sorted stream. For those, run `bai-build` on the output instead. The index is byte-identical to `bai-build`. |

The thread/memory knobs (except `--compression`) apply to the sorted-merge path
only; `--concat` is single-threaded BGZF passthrough.

`--concat` has priority over `--order`. With both options, the tool
concatenates the inputs and does not show an error. The `--order` value is
still checked first: an invalid order, or `--order tag` without `--tag`, is an
error.

On stderr, the tool prints a summary line, then a runtime report. The report
shows the wall time, the peak RSS, the CPU time, and the time of each pipeline
stage. It ends with a `verdict` line that names the probable bottleneck. The
`--concat` and disjoint-chain passthrough paths have no stages: their report
says `pipeline: BGZF block passthrough` and `verdict: I/O-COPY bound`.

### bai-query — region query

Extract records overlapping a genomic region using a BAI index.
Expects a `.bai` file alongside the BAM (auto-detected as `<bam>.bai`).

```sh
pbsamoa bai-query sorted.bam chr1:1000-2000
```

Arguments: BAM path, region (`ref`, `ref:start`, or `ref:start-end`, 1-based
inclusive; split on the rightmost colon, as in `dump --region`).

The output is SAM text on stdout: the header, then the records. If
`<bam>.bai` does not exist, the tool stops with an error. A reference that is
not in the header is an error too, but the tool has already written the header
to stdout by then.

### bai-build — build BAI

Build a `.bai` index from a coordinate-sorted BAM file.

```sh
pbsamoa bai-build sorted.bam
# creates sorted.bam.bai
pbsamoa bai-build --threads 4 sorted.bam
```

Options:

| Flag | Description |
| ---- | ----------- |
| `--threads N` | BGZF decompression workers; `0` = auto, `min(hw, 8)` (default: `0`) |

On stderr the tool prints `Index written to <path>`.

### zmi-build — build ZMW index

Copy a BAM file and build a `.zmi` ZMW index alongside it.

```sh
pbsamoa zmi-build input.bam output.bam
# creates output.bam and output.bam.zmi
```

On stderr the tool prints `ZMI index written to <path>`.

Both ZMI builders capture ordered barcode IDs from the BAM `bc` tag in ZMI 1.1.
The library's `ZmwIndex::FindBarcodes()` can retrieve matching record offsets;
`zmi-query` continues to query by ZMW.

### zmi-index — build `.zmi` for an existing BAM

Build a `.zmi` sidecar from an existing BAM without rewriting the BAM.
Uses parallel BGZF decompression. `--threads` defaults to
`hardware_concurrency()`.

```sh
pbsamoa zmi-index input.bam
# creates input.bam.zmi

pbsamoa zmi-index --threads 8 --quiet input.bam
```

Options:

| Flag | Description |
| ---- | ----------- |
| `--threads N` | Worker threads; `0` = auto, all hardware threads (default: `0`). A negative value is an error |
| `--quiet` | Do not write the summary line `ZMI index written to <path> (N records)` to stderr |

### zmi-query — query by ZMW

Extract the records of one ZMW hole number.

```sh
pbsamoa zmi-query movie.bam 12345
```

The output is SAM text on stdout: the header, then the records. The tool needs
a ZMW index (`.zmi` or `.pbi`) next to the BAM. Without one, it stops with the
error `No ZMW index found for: <bam>`.

### chunk — partitioned reading

Read a partition of a BAM file, divided by unique ZMWs. Records of a ZMW always
stay together (works for subreads and HiFi), and running all `TOTAL` chunks
visits every record exactly once.

```sh
pbsamoa chunk movie.bam 1 4   # chunk 1 of 4 (1-based)
pbsamoa chunk movie.bam 4 4   # chunk 4 of 4 (last)
```

`CHUNK` must be in `[1, TOTAL]` and `TOTAL` at least 1; other values are an
error (status `1`).

By default each chunk is one contiguous ZMW range. Use `--mode scatter` for a
chaotic-but-deterministic partition that samples ZMWs across the whole file —
useful when a single chunk should be a representative sample of the dataset:

```sh
pbsamoa chunk movie.bam 1 4 --mode scatter            # spread across the file
pbsamoa chunk movie.bam 1 4 --mode scatter --tile 64  # 64 ZMWs/seek (fewer seeks)
pbsamoa chunk movie.bam 1 4 --mode scatter --seed 7   # reproducible shuffle seed
```

Scatter assigns tiles of up to `M` consecutive ZMWs (`--tile M`, default 100) to
chunks via a seeded balanced shuffle (`--seed S`, default 42). The same
`(file, CHUNK, TOTAL, M, S)` always yields the same chunk on every platform.
Larger `M` reads more ZMWs per seek (less disk seeking, coarser sampling).
Both modes require a ZMW index (`.zmi`/`.pbi`); see `zmi-build`.

The output is SAM text on stdout: the header, then the records of the chunk.
`--tile` and `--seed` without `--mode scatter` are an error.

Options:

| Flag | Description |
| ---- | ----------- |
| `--mode MODE` | `contiguous` \| `scatter` (default: `contiguous`) |
| `--tile M` | Scatter only: up to `M` consecutive ZMWs per seek, `>= 1` (default: `100`) |
| `--seed S` | Scatter only: unsigned 64-bit shuffle seed (default: `42`) |

### bench — benchmarks

Run throughput benchmarks on a BAM file. Executes multiple benchmarks
in sequence: `sequential_read`, `batch_read`, `pipeline_raw_read`,
`record_reader` (serial), `record_reader` (parallel), `write`, and
`region_query`.

```sh
pbsamoa bench input.bam
```

Options:

| Flag | Description |
| ---- | ----------- |
| `--bgzf-threads N` | BGZF decompression worker threads; `-1` = auto, `min(hw, 8)` (default: `-1`). An explicit count is capped at 10 |
| `--decode-threads N` | BamRecord decode worker threads; `0` = serial (default: 4) |

With `--bgzf-threads 0`, the tool skips `pipeline_raw_read`. With
`--bgzf-threads 0 --decode-threads 0`, it also skips the parallel
`record_reader`. Without `<input>.bai`, it skips `region_query`. The `write`
benchmark writes a temporary BAM in the system temp directory and removes it.
The tool writes the results to stdout.
