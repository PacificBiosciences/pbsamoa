/* crc_base_aliases.c of ISA-L 2.31.1 (BSD-3-Clause, Copyright(c) 2011-2017 Intel Corporation),
 * with one change for pbsamoa: crc32_gzip_refl selects the PCLMUL kernel at run time. */

#include "crc.h"
#include "crc64.h"
#include <stdint.h>

unsigned int
crc32_iscsi(unsigned char *buffer, int len, unsigned int crc_init)
{
        return crc32_iscsi_base(buffer, len, crc_init);
}

uint16_t
crc16_t10dif(uint16_t seed, const unsigned char *buf, uint64_t len)
{
        return crc16_t10dif_base(seed, (uint8_t *) buf, len);
}

uint16_t
crc16_t10dif_copy(uint16_t seed, uint8_t *dst, uint8_t *src, uint64_t len)
{
        return crc16_t10dif_copy_base(seed, dst, src, len);
}

uint32_t
crc32_ieee(uint32_t seed, const unsigned char *buf, uint64_t len)
{
        return crc32_ieee_base(seed, (uint8_t *) buf, len);
}

uint32_t
crc32_gzip_refl_by8(uint32_t seed, const unsigned char *buf, uint64_t len);

/* The gzip CRC of igzip inflate: the PCLMUL kernel of the upstream dispatcher on CPUs with
 * SSE4.1 and PCLMUL, else the portable code (byte by byte: with it, inflate of a FASTQ stream is
 * slower than zlib) */
uint32_t
crc32_gzip_refl(uint32_t seed, const unsigned char *buf, uint64_t len)
{
        if (__builtin_cpu_supports("pclmul") && __builtin_cpu_supports("sse4.1"))
                return crc32_gzip_refl_by8(seed, buf, len);
        return crc32_gzip_refl_base(seed, (uint8_t *) buf, len);
}

uint64_t
crc64_ecma_refl(uint64_t seed, const uint8_t *buf, uint64_t len)
{
        return crc64_ecma_refl_base(seed, buf, len);
}

uint64_t
crc64_ecma_norm(uint64_t seed, const uint8_t *buf, uint64_t len)
{
        return crc64_ecma_norm_base(seed, buf, len);
}

uint64_t
crc64_iso_refl(uint64_t seed, const uint8_t *buf, uint64_t len)
{
        return crc64_iso_refl_base(seed, buf, len);
}

uint64_t
crc64_iso_norm(uint64_t seed, const uint8_t *buf, uint64_t len)
{
        return crc64_iso_norm_base(seed, buf, len);
}

uint64_t
crc64_jones_refl(uint64_t seed, const uint8_t *buf, uint64_t len)
{
        return crc64_jones_refl_base(seed, buf, len);
}

uint64_t
crc64_jones_norm(uint64_t seed, const uint8_t *buf, uint64_t len)
{
        return crc64_jones_norm_base(seed, buf, len);
}

uint64_t
crc64_rocksoft_refl(uint64_t seed, const uint8_t *buf, uint64_t len)
{
        return crc64_rocksoft_refl_base(seed, buf, len);
}

uint64_t
crc64_rocksoft_norm(uint64_t seed, const uint8_t *buf, uint64_t len)
{
        return crc64_rocksoft_norm_base(seed, buf, len);
}
