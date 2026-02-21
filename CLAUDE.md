# CLAUDE.md

pbsamoa — C++23 SAM/BAM/BAI library. Single dep: libdeflate (Meson wrap).

## Build

Standard meson (see `cpp.md`). Additional:
```sh
# clang-tidy
rg --files -g 'src/**/*.cpp' -g 'src/**/*.hpp' -g 'include/**/*.hpp' -g 'tests/unit/**/*.cpp' \
  | xargs -n1 -P8 clang-tidy -p build --extra-arg-before=--config=/opt/homebrew/etc/clang/arm64-apple-darwin25.cfg
```

## Architecture (four layers, bottom-up)

**Compression** — `VirtualOffset`, `BgzfReader`, `BgzfWriter`, `BgzfPipeline` (3-stage parallel decompression via libdeflate).

**Record** — `RawRecord` (owning, copies raw BAM bytes, inline decode-on-demand), `RawRecordBatch` (owns decompressed buffers, indexed access), `BamRecord` (structured fields, zero-cost mutation, serialize once in writer), `CigarOp`, `TagMap`.

**API** — `BamRawReader` (zero-copy iteration, sync/parallel BGZF), `BamRecordReader` (wraps BamRawReader, pre-decodes to `BamRecord` in background thread), `SamReader`, `BamWriter`, `SamWriter`, `BaiIndex`. Both readers: range-based + batch interfaces, memory-budget sizing, region queries via BAI.

**ZMW Index** — `ZmwIndex` (`.zmi`/`.pbi`, point/batch queries, chunking by `(rgId, zmw)`), `ZmwIdentity`, `ZmiWriter`, `ZmiBamWriter` (BAM+ZMI coordinated output), `BamZmwReader` (groups by ZMW). Chunking via `BamRawReaderConfig::RecordLimit`.

## Design

- `RawRecord`: decode logic in header for inlining
- `BamRecord`: structured data (not BAM layout) — serialize once in writer
- `ToOwned()` tag filtering: `DropTags`/`KeepTags` (skip PacBio `ip`/`pw`)
- Batches by `ByteLimit`, not record count
- Writer accepts `BamRecord`, `RawRecord`, `span<const byte>`

## Testing

Cram tests use `samtools` as reference. Patterns:
```
  $ pbsamoa-dump "${TESTDIR}"/../data/small.bam
  @HD\tVN:1.6\tSO:coordinate (esc)
  *r001*99*ref* (glob)
```
- `(glob)` wildcard, `(re)` regex, `(esc)` escapes, `[1]` expected non-zero exit
- Fixtures: `"${TESTDIR}"/../data/`; tool binaries via env vars in `tests/meson.build`

## Reference

`reference/spec/SAMv1.md`, `reference/spec/SAMtags.md`, `reference/htslib/`, `docs/architecture.md`, `docs/api.md`, `docs/cli.md`
