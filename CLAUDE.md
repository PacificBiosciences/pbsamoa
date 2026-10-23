# CLAUDE.md

pbsamoa — C++23 SAM/BAM/CRAM/BAI library. Core deps (Meson wraps): libdeflate + pbcopper (pbcopper also powers the CLI tools via CLIv2; it needs Boost from the system). Optional deps: bzip2 and liblzma (system, for CRAM bzip2/lzma blocks) and ISA-L (x86-64 with NASM, option `isal`). The CRAM codecs come from the vendored `third-party/htscodecs`. `AGENTS.md` points non-Claude agents to this file.

## Build & Test

Standard meson:
```sh
meson setup build && meson compile -C build && meson test -C build
```

clang-tidy across project sources (config: `.clang-tidy`):
```sh
rg --files -g 'src/**/*.cpp' -g 'src/**/*.hpp' -g 'include/**/*.hpp' -g 'tests/unit/**/*.cpp' \
  | xargs -n1 -P8 clang-tidy -p build
```
With Homebrew LLVM on macOS, add `--extra-arg-before=--config=$(brew --prefix llvm)/etc/clang/<target>.cfg`.

Cram tests get the tool path from env vars set in `tests/cram/meson.build` (`PBSAMOA`, `PATH`) — run via `meson test`, not direct invocation.

## Architecture (`include/pbsamoa/{core,io,index,cram}`; `docs/architecture.md` groups them as Compression, Record, CRAM, Indexes, API)

**core/ — Compression** — `VirtualOffset`, `BgzfReader` (sync or parallel via `numWorkers`, 3-stage pipeline internally), `BgzfWriter`, `BgzfWriterPool`, `BgzfBlockEncoder` + `BgzfBlockWriter` (callers compress on their own threads; the writer only appends blocks).

**core/ — Record** — `RawRecordView` (non-owning, inline decode-on-demand over borrowed bytes), `RawRecord` (owns a copy of the record bytes; borrow via `View()`), `RawRecordBatch` (owns decompressed buffers; `View(i)` for copy-free access), `BamRecord` (structured fields, zero-cost mutation, serialize once in writer), `CigarOp`, `TagMap`, `SamHeader`, `GenomicInterval`, `SequenceView` (4-bit base codes, `Sequence.hpp`), `Basemods` (MM/ML tags). Clipping: `CigarClipping` (query or reference coordinates), `TagClipper` + `TagClipStrategy` subclasses (per-tag strategies, `TagClipping.hpp`), `LosslessClipping` (keeps the clipped flanks so `undo` can restore them).

**io/ — BAM/SAM/CRAM API** — `BamRawReader` (sync/parallel BGZF, range + batch + whitelist interfaces, region queries via BAI; `Records()` yields owning `RawRecord`, `ReadBatch()` + `View(i)` is the copy-free path), `BamRecordReader` (wraps BamRawReader, pre-decodes to `BamRecord` via background producer thread + optional thread pool), `SamReader`, `BamWriter`, `SamWriter`, `BamFile` (path + header + `.bai`), `BamCollection` (ordered inputs with a merged header). `BamSort.hpp` (`SortConfig` + `SortBam`, external RAM-bounded sort), `BamMerge.hpp` (`MergeConfig` + `MergeBam`, k-way merge; implemented in `src/BamSort.cpp`) + `MergeReadAhead` (per-source read-ahead; output is independent of thread count). `CramReader`, `CramWriter` (reference optional; default output is reference-free). ZMW: `BamZmwReader` (groups by ZMW), `ZmiWriter`, `ZmiBamWriter` (BAM+ZMI coordinated output). Memory-budget sizing via `ByteLimit`.

**index/** — `BaiIndex`, `BaiStreamBuilder` (BAI on the fly from writer callbacks), `CraiIndex` (CRAI interop), `ZmwIndex` (`.zmi`/`.pbi`, point/batch queries, chunking by `(rgId, zmw)`), `ZmwIdentity`, `ZmwWhitelist`. Chunking via `BamRawReaderConfig::ChunkNum`/`TotalChunks`.

**cram/** — `CramCodec` (data-series codecs), `CramCompression.hpp` (`CompressionRegistry`: block methods incl. FQZ slices), `CramStructs`, `CramReference` (the reference in windows, for the writer and the reader). Round-trips with `BamRecord`.

## CLI Tools (`src/tools/<name>/`, one `pbsamoa` executable from `src/tools/pbsamoa.cpp`)

`pbsamoa` (umbrella) · `dump` · `convert` · `sort` · `merge` · `subsample` · `depth` · `chunk` · `undo` · `bai-build` · `bai-query` · `zmi-build` · `zmi-index` · `zmi-query` · `bench`. Shared helpers: `CliUtils.hpp`, `ParseUtils.hpp`, `SamOutput.hpp`, `MetricUtils.hpp`, `StatusStream.{hpp,cpp}`. `merge` auto-detects sorted (k-way merge) vs unsorted (`--concat`, BGZF block passthrough).

## Design

- `RawRecord`: decode logic in header for inlining
- `BamRecord`: structured data (not BAM layout) — serialize once in writer
- `ToOwned()` tag filtering: `DropTags`/`KeepTags` (skip PacBio `ip`/`pw`)
- Batches by `ByteLimit`, not record count
- Writer accepts `BamRecord`, `RawRecord`, `span<const byte>`

## Cram test syntax

```
  $ "${PBSAMOA}" dump "${TESTDIR}"/../data/spec_example.bam
  @HD\tVN:1.6\tSO:coordinate (esc)
  *r001*99*ref* (glob)
```
`(glob)` wildcard · `(re)` regex · `(esc)` escapes · `[N]` expected exit status N. Fixtures: `tests/data/`. Reference oracle: `samtools`.

## Reference

`reference/spec/SAMv1.md`, `reference/spec/SAMtags.md`, `reference/spec/CRAMv3.md`, `reference/htslib/` (ignored; clone htslib there for a local copy), `docs/architecture.md`, `docs/api.md`, `docs/cli.md`
