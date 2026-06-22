# CLAUDE.md

pbsamoa — C++23 SAM/BAM/CRAM/BAI library. Single dep: libdeflate (Meson wrap). `AGENTS.md` is the non-Claude mirror; keep both in sync.

## Build & Test

Standard meson (see `cpp.md`):
```sh
meson setup build && meson compile -C build && meson test -C build
```

clang-tidy across project sources:
```sh
rg --files -g 'src/**/*.cpp' -g 'src/**/*.hpp' -g 'include/**/*.hpp' -g 'tests/unit/**/*.cpp' \
  | xargs -n1 -P8 clang-tidy -p build --extra-arg-before=--config=/opt/homebrew/etc/clang/arm64-apple-darwin25.cfg
```

Cram tests resolve tool paths via env vars in `tests/meson.build` — run via `meson test`, not direct invocation.

## Architecture (five layers, bottom-up)

**Compression** — `VirtualOffset`, `BgzfReader` (sync or parallel via `numWorkers`, 3-stage pipeline internally), `BgzfWriter`.

**Record** — `RawRecord` (owning, copies raw BAM bytes, inline decode-on-demand), `RawRecordBatch` (owns decompressed buffers, indexed access), `BamRecord` (structured fields, zero-cost mutation, serialize once in writer), `CigarOp`, `TagMap`.

**BAM/SAM API** — `BamRawReader` (zero-copy iteration, sync/parallel BGZF, range + batch + whitelist interfaces, region queries via BAI), `BamRecordReader` (wraps BamRawReader, pre-decodes to `BamRecord` via background producer thread + optional thread pool), `SamReader`, `BamWriter`, `SamWriter`, `BaiIndex`. Memory-budget sizing via `ByteLimit`.

**CRAM v3** — `CramReader`, `CramWriter`, `CraiIndex` (CRAI interop), `CramCodec`/`CramCompression` (configurable block + data-series codecs incl. FQZ slices), `CramStructs`, `CramMd5`. Reference-based; round-trips with `BamRecord`.

**ZMW Index** — `ZmwIndex` (`.zmi`/`.pbi`, point/batch queries, chunking by `(rgId, zmw)`), `ZmwIdentity`, `ZmwWhitelist`, `ZmiWriter`, `ZmiBamWriter` (BAM+ZMI coordinated output), `BamZmwReader` (groups by ZMW). Chunking via `BamRawReaderConfig::ChunkNum`/`TotalChunks`.

## CLI Tools (`src/tools/<name>/`)

`pbsamoa` (umbrella) · `dump` · `convert` · `sort` · `merge` · `chunk` · `bai-build` · `bai-query` · `zmi-build` · `zmi-index` · `zmi-query` · `bench`. Shared helpers: `CliUtils.hpp`, `ParseUtils.hpp`, `SamOutput.hpp`, `MetricUtils.hpp`. `merge` auto-detects sorted (k-way merge) vs unsorted (`--concat`, BGZF block passthrough).

## Design

- `RawRecord`: decode logic in header for inlining
- `BamRecord`: structured data (not BAM layout) — serialize once in writer
- `ToOwned()` tag filtering: `DropTags`/`KeepTags` (skip PacBio `ip`/`pw`)
- Batches by `ByteLimit`, not record count
- Writer accepts `BamRecord`, `RawRecord`, `span<const byte>`

## Cram test syntax

```
  $ pbsamoa-dump "${TESTDIR}"/../data/small.bam
  @HD\tVN:1.6\tSO:coordinate (esc)
  *r001*99*ref* (glob)
```
`(glob)` wildcard · `(re)` regex · `(esc)` escapes · `[1]` expected non-zero exit. Fixtures: `tests/data/`. Reference oracle: `samtools`.

## Reference

`reference/spec/SAMv1.md`, `reference/spec/SAMtags.md`, `reference/htslib/`, `docs/architecture.md`, `docs/api.md`, `docs/cli.md`
