# CRAM Format Specification (Version 3.0/3.1)

> Converted from the [samtools hts-specs CRAMv3.tex](https://github.com/samtools/hts-specs).
> License: Apache 2.0

## 1. Overview

CRAM has the following major objectives:

1. Significantly better lossless compression than BAM
2. Full compatibility with BAM
3. Effortless transition to CRAM from using BAM files
4. Support for controlled loss of BAM data

Data in CRAM is stored either as CRAM records or using general purpose compressors
(gzip, bzip2, lzma, rANS). CRAM records are compressed using a number of different
encoding strategies. For example, bases are reference compressed by encoding base
differences rather than storing the bases themselves.

## 2. Data Types

### 2.1 Logical Data Types

| Type    | Description              |
|---------|--------------------------|
| Byte    | Signed byte (8 bits)     |
| Integer | Signed 32-bit integer    |
| Long    | Signed 64-bit integer    |
| Array   | Array of any logical type: `array<type>` |

### 2.2 Storage Data Types (Byte Stream)

CRAM uses **little endianness** for bytes.

| Type          | Code   | Description                                      |
|---------------|--------|--------------------------------------------------|
| Boolean       | bool   | 1 byte: 0x0 = false, 0x1 = true                 |
| Integer       | int32  | Signed 32-bit, 4 bytes little-endian             |
| Long          | int64  | Signed 64-bit, 8 bytes little-endian             |
| ITF-8         | itf8   | Variable-length integer (1-5 bytes, like UTF-8)  |
| LTF-8         | ltf8   | Variable-length long (1-9 bytes)                 |
| Array         | array  | ITF-8 length prefix followed by elements         |
| Encoding      | —      | ITF-8 codec_id + ITF-8 param_length + params     |
| Map           | —      | ITF-8 size_bytes + ITF-8 count + key/value pairs |
| String        | —      | UTF-8 byte arrays                                |

### 2.3 ITF-8 Encoding

Prefix bits determine number of following bytes:
- `0xxxxxxx` — 1 byte total (7 data bits)
- `10xxxxxx` + 1 byte — 2 bytes total (14 data bits)
- `110xxxxx` + 2 bytes — 3 bytes total (21 data bits)
- `1110xxxx` + 3 bytes — 4 bytes total (28 data bits)
- `11110xxx` + 4 bytes — 5 bytes total (32 data bits, last byte uses 4 bits)

### 2.4 Bit Stream

The CORE block supports bit-based encoding methods. Bits are written most
significant bit first. The last byte is left-shifted to fill remaining bits.

## 3. Encodings

| Codec          | ID | Parameters                               | Description                      |
|----------------|----|------------------------------------------|----------------------------------|
| NULL           | 0  | none                                     | Series not preserved             |
| EXTERNAL       | 1  | int block_content_id                     | Data in external block           |
| _(GOLOMB)_     | 2  | int offset, int M                        | Deprecated                       |
| HUFFMAN        | 3  | array\<int\>, array\<int\>               | Huffman coding                   |
| BYTE_ARRAY_LEN | 4  | encoding\<int\> len, encoding\<byte\> bytes | Byte arrays with length       |
| BYTE_ARRAY_STOP| 5  | byte stop, int external_block_id         | Byte arrays with stop value      |
| BETA           | 6  | int offset, int num_bits                 | Binary coding                    |
| SUBEXP         | 7  | int offset, int K                        | Subexponential coding            |
| _(GOLOMB_RICE)_| 8  | int offset, int log2m                    | Deprecated                       |
| GAMMA          | 9  | int offset                               | Elias gamma coding               |

## 4. Checksums

- **CRC32**: Polynomial 0x04C11DB7 (ITU-T V.42). Written as int32.
- **CRC32 sum**: Sum of individual CRC32 values modulo 2^32.

## 5. File Structure

```
┌─────────────────┬──────────────────────┬────────────────┬─────┬────────────────┬──────────────────────┐
│ File Definition │ CRAM Header Container│ Data Container │ ... │ Data Container │ CRAM EOF Container   │
│ (26 bytes)      │                      │                │     │                │                      │
└─────────────────┴──────────────────────┴────────────────┴─────┴────────────────┴──────────────────────┘
```

### 5.1 File Definition (26 bytes)

| Data type      | Name                 | Value                                    |
|----------------|----------------------|------------------------------------------|
| byte[4]        | magic                | "CRAM" (0x43 0x52 0x41 0x4d)            |
| unsigned byte  | major version        | 3                                        |
| unsigned byte  | minor version        | 0 or 1                                   |
| byte[20]       | file id              | File identifier (e.g. filename or SHA1)  |

### 5.2 Container Header

| Data type     | Name               | Description                                        |
|---------------|--------------------|----------------------------------------------------|
| int32         | length             | Total byte length of blocks + padding               |
| itf8          | ref_seq_id         | Reference ID, -1=unmapped, -2=multi-ref            |
| itf8          | start_pos          | Alignment start position                           |
| itf8          | alignment_span     | Length of alignment                                |
| itf8          | num_records        | Number of records in container                     |
| ltf8          | record_counter     | 0-based sequential record index                    |
| ltf8          | bases              | Number of read bases                               |
| itf8          | num_blocks         | Total number of blocks                             |
| array\<itf8\> | landmarks          | Byte offsets of slices from end of container header |
| int32         | crc32              | CRC32 of preceding bytes                           |

### 5.3 Block Structure

| Data type | Name              | Description                                        |
|-----------|-------------------|----------------------------------------------------|
| byte      | method            | Compression: 0=raw, 1=gzip, 2=bzip2, 3=lzma, 4=rans4x8, 5=rans4x16, 6=arith, 7=fqzcomp, 8=tok |
| byte      | content_type_id   | 0=FILE_HEADER, 1=COMPRESSION_HEADER, 2=SLICE_HEADER, 4=EXTERNAL_DATA, 5=CORE_DATA |
| itf8      | content_id        | Block content identifier                           |
| itf8      | compressed_size   | Size after compression                             |
| itf8      | raw_size          | Size before compression                            |
| byte[]    | data              | Block data (compressed)                            |
| byte[4]   | crc32             | CRC32 of all preceding bytes in block              |

### 5.4 Compression Header Block

Three parts:

1. **Preservation Map** (byte[2] keys):
   - `RN`: Read names included (bool, default true)
   - `AP`: AP data series is delta-encoded (bool, default true)
   - `RR`: Reference required (bool, default true)
   - `SM`: Substitution matrix (byte[5])
   - `TD`: Tag IDs dictionary (array\<byte\>)

2. **Data Series Encoding Map** (byte[2] keys → encoding):
   - `BF`: BAM bit flags, `CF`: CRAM bit flags
   - `RI`: Reference ID, `RL`: Read length, `AP`: Alignment position
   - `RG`: Read group, `RN`: Read names
   - `MF`: Mate flags, `NS`: Next mate ref, `NP`: Next mate pos, `TS`: Template size
   - `NF`: Distance to next fragment, `TL`: Tag IDs index
   - `FN`: Number of read features, `FC`: Feature codes, `FP`: Feature positions
   - `DL`: Deletion lengths, `BB`: Base stretches, `QQ`: Quality stretches
   - `BS`: Base substitution codes, `IN`: Insertions
   - `RS`: Reference skip, `PD`: Padding, `HC`: Hard clip, `SC`: Soft clip
   - `MQ`: Mapping qualities, `BA`: Bases, `QS`: Quality scores

3. **Tag Encoding Map**: Keys are ITF-8 of `(char1<<16)+(char2<<8)+type`

### 5.5 Slice Header Block

| Data type | Name                        | Description                               |
|-----------|-----------------------------|-------------------------------------------|
| itf8      | ref_seq_id                  | Reference ID (-1=unmapped, -2=multi-ref)  |
| itf8      | alignment_start             | Alignment start position                  |
| itf8      | alignment_span              | Length of alignment                       |
| itf8      | num_records                 | Number of records in slice                |
| ltf8      | record_counter              | 0-based sequential record index           |
| itf8      | num_blocks                  | Number of blocks in slice                 |
| itf8[]    | block_content_ids           | Content IDs of blocks in slice            |
| itf8      | embedded_ref_block_id       | Block content ID for embedded ref, or -1  |
| byte[16]  | ref_md5                     | MD5 of reference within slice boundaries  |
| byte[]    | optional_tags               | BAM-encoded auxiliary tags                |

## 6. EOF Container (38 bytes)

```
0f 00 00 00 ff ff ff ff 0f e0 45 4f 46 00 00 00 00 01 00 05 bd d9 4f
00 01 00 06 06 01 00 01 00 01 00 ee 63 01 4b
```

## 7. Record Structure

### 7.1 Common Fields (all records)

1. **BF** — BAM bit flags (int)
2. **CF** — CRAM bit flags (int): 0x1=qual_as_array, 0x2=detached, 0x4=has_mate_downstream, 0x8=decode_seq_as_star
3. Positional data (RI, RL, AP, RG)
4. Read names (if RN preserved or CF & 0x2)
5. Mate data (if CF & 0x2: MF, NS, NP, TS; if CF & 0x4: NF)
6. Tag data (TL, then tag values per TD dictionary)
7. Mapped read features OR unmapped bases/qualities

### 7.2 Read Feature Codes

| Code | Name              | Data series          |
|------|-------------------|----------------------|
| B    | Base substitution | FC, FP, BS           |
| X    | Substitution      | FC, FP, BS           |
| I    | Insertion         | FC, FP, IN           |
| D    | Deletion          | FC, FP, DL           |
| i    | Single insert     | FC, FP, BA           |
| b    | Stretch of bases  | FC, FP, BB           |
| q    | Quality score     | FC, FP, QS           |
| Q    | Stretch of quals  | FC, FP, QQ           |
| S    | Soft clip         | FC, FP, SC           |
| N    | Reference skip    | FC, FP, RS           |
| P    | Padding           | FC, FP, PD           |
| H    | Hard clip         | FC, FP, HC           |

### 7.3 Mapped Read Decode

For mapped reads (BAM flag 0x4 not set):
1. Decode read features (FN count, then FC/FP + feature-specific data series)
2. Decode mapping quality (MQ)
3. Decode quality scores (if CF & 0x1, decode RL quality scores from QS)

### 7.4 Unmapped Read Decode

For unmapped reads (BAM flag 0x4 set):
1. Decode bases: RL values from BA data series
2. Decode quality scores: if CF & 0x1, decode RL values from QS data series

## 8. Compression Methods (Block Level)

| ID | Method     | CRAM Version | Description                        |
|----|------------|--------------|------------------------------------|
| 0  | raw        | all          | No compression                     |
| 1  | gzip       | all          | RFC 1952                           |
| 2  | bzip2      | 2.0+         | BWT-based compression              |
| 3  | lzma       | 3.0+         | Lempel-Ziv-Markov chain            |
| 4  | rans4x8    | 3.0+         | rANS with 4 interleaved streams    |
| 5  | rans4x16   | 3.1          | rANS with 16-bit renormalization   |
| 6  | arith      | 3.1          | Adaptive arithmetic coder          |
| 7  | fqzcomp    | 3.1          | Quality score compressor           |
| 8  | tok        | 3.1          | Name tokeniser                     |

## 9. Version Differences

- **3.0**: Introduced LZMA and rANS4x8 codecs, header/data CRC32 checksums, improvements for unsorted data.
- **3.1**: Added rANS4x16, adaptive arithmetic, fqzcomp, and name tokeniser codecs. Otherwise identical to 3.0.

Tools that don't use 3.1 codecs should write version 3.0 for maximum compatibility.

## 10. Reference

- Full specification: https://samtools.github.io/hts-specs/CRAMv3.pdf
- Source: https://github.com/samtools/hts-specs
