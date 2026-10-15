#!/usr/bin/env python3
"""Write BAM.pbi (PacBio index, BasicData only) for an unaligned PacBio BAM.

Usage: make_pbi.py BAM

Each row holds rgId, qStart = 0, qEnd = read length, holeNumber, readQual,
ctxtFlag = 0, and the virtual offset of the record, as pbindex writes them for
unaligned HiFi reads.
"""

import bisect
import struct
import sys
import zlib

BGZF_EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")


def read_bgzf_blocks(path):
    """Yield (compressed offset, decompressed bytes) for each BGZF block."""
    with open(path, "rb") as handle:
        data = handle.read()
    pos = 0
    while pos < len(data):
        block_size = struct.unpack_from("<H", data, pos + 16)[0] + 1
        payload = data[pos + 18 : pos + block_size - 8]
        yield pos, zlib.decompress(payload, -15)
        pos += block_size


def records(path):
    """Yield (virtual offset, record body) for each BAM record."""
    stream = bytearray()
    starts = []  # (offset in stream, compressed block offset)
    for coffset, block in read_bgzf_blocks(path):
        starts.append((len(stream), coffset))
        stream += block

    block_starts = [start for start, _ in starts]

    def virtual_offset(pos):
        block = bisect.bisect_right(block_starts, pos) - 1
        start, coffset = starts[block]
        return (coffset << 16) | (pos - start)

    text_length = struct.unpack_from("<i", stream, 4)[0]
    pos = 8 + text_length
    num_refs = struct.unpack_from("<i", stream, pos)[0]
    pos += 4
    for _ in range(num_refs):
        name_length = struct.unpack_from("<i", stream, pos)[0]
        pos += 4 + name_length + 4
    while pos < len(stream):
        block_size = struct.unpack_from("<i", stream, pos)[0]
        yield virtual_offset(pos), bytes(stream[pos + 4 : pos + 4 + block_size])
        pos += 4 + block_size


def find_z_tag(aux, key):
    """Return the value of a Z tag, or None. Handles the tag types of this fixture."""
    sizes = {"A": 1, "c": 1, "C": 1, "s": 2, "S": 2, "i": 4, "I": 4, "f": 4}
    pos = 0
    while pos < len(aux):
        tag, kind = aux[pos : pos + 2].decode(), chr(aux[pos + 2])
        pos += 3
        if kind in "ZH":
            end = aux.index(0, pos)
            if tag == key:
                return aux[pos:end].decode()
            pos = end + 1
        elif kind == "B":
            sub, count = chr(aux[pos]), struct.unpack_from("<i", aux, pos + 1)[0]
            pos += 5 + sizes[sub] * count
        else:
            pos += sizes[kind]
    return None


def bgzf_compress(data):
    out = bytearray()
    for start in range(0, len(data), 0xFF00):
        chunk = data[start : start + 0xFF00]
        compressor = zlib.compressobj(6, zlib.DEFLATED, -15)
        deflated = compressor.compress(chunk) + compressor.flush()
        header = struct.pack("<4BI2BH2BHH", 0x1F, 0x8B, 8, 4, 0, 0, 0xFF, 6, 66, 67, 2,
                             len(deflated) + 25)
        out += header + deflated + struct.pack("<II", zlib.crc32(chunk), len(chunk))
    return bytes(out) + BGZF_EOF


def main():
    bam = sys.argv[1]
    rows = []
    for voffset, body in records(bam):
        name_length = body[8]
        cigar_ops = struct.unpack_from("<H", body, 12)[0]
        seq_length = struct.unpack_from("<i", body, 16)[0]
        name = body[32 : 32 + name_length - 1].decode()
        aux = body[32 + name_length + 4 * cigar_ops + (seq_length + 1) // 2 + seq_length :]
        rg = find_z_tag(aux, "RG")
        rg_id = struct.unpack("<i", struct.pack("<I", int(rg.split("/")[0], 16)))[0]
        rows.append((rg_id, seq_length, int(name.split("/")[1]), voffset))

    n = len(rows)
    out = bytearray(b"PBI\x01" + struct.pack("<IHI", 0x030001, 0, n) + bytes(18))
    out += struct.pack(f"<{n}i", *(r[0] for r in rows))
    out += struct.pack(f"<{n}i", *([0] * n))
    out += struct.pack(f"<{n}i", *(r[1] for r in rows))
    out += struct.pack(f"<{n}i", *(r[2] for r in rows))
    out += struct.pack(f"<{n}f", *([0.99] * n))
    out += bytes(n)
    out += struct.pack(f"<{n}q", *(r[3] for r in rows))
    with open(bam + ".pbi", "wb") as handle:
        handle.write(bgzf_compress(bytes(out)))


if __name__ == "__main__":
    main()
