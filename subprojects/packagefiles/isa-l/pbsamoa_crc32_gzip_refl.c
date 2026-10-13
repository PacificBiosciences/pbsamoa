/* The one CRC that igzip calls, for pbsamoa (from crc_base_aliases.c of ISA-L, BSD-3-Clause,
 * Copyright(c) 2011-2017 Intel Corporation). The upstream build gets it from crc_multibinary.asm,
 * which needs every CRC kernel of ISA-L; pbsamoa builds only igzip. */

#include "crc.h"
#include <stdint.h>

uint32_t
crc32_gzip_refl_by8(uint32_t seed, const unsigned char *buf, uint64_t len);

/* The gzip CRC of igzip deflate and inflate: the PCLMUL kernel of the upstream dispatcher on CPUs
 * with SSE4.1 and PCLMUL, else the portable code (byte by byte: with it, inflate of a FASTQ stream
 * is slower than zlib) */
uint32_t
crc32_gzip_refl(uint32_t seed, const unsigned char *buf, uint64_t len)
{
        if (__builtin_cpu_supports("pclmul") && __builtin_cpu_supports("sse4.1"))
                return crc32_gzip_refl_by8(seed, buf, len);
        return crc32_gzip_refl_base(seed, (uint8_t *) buf, len);
}
