#include <stdint.h>
#include <string.h>
#include <immintrin.h>

uint64_t hts_pack_avx2_count_symbols(const uint8_t *data, uint64_t len, int *present) {
    uint8_t symbols[16];
    int n = 0;
    uint64_t i;
    for (i = 0; i + 32 <= len; i += 32) {
        __m256i bytes = _mm256_loadu_si256((const __m256i *)(data + i));
        __m256i matched = _mm256_setzero_si256();
        for (int s = 0; s < n; ++s)
            matched = _mm256_or_si256(matched, _mm256_cmpeq_epi8(bytes, _mm256_set1_epi8(symbols[s])));
        uint32_t missing = ~(uint32_t)_mm256_movemask_epi8(matched);
        while (missing) {
            uint8_t value = data[i + __builtin_ctz(missing)];
            if (n == 16)
                return i;
            present[value] = 1;
            symbols[n++] = value;
            missing &= ~(uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(bytes, _mm256_set1_epi8(value)));
        }
    }
    return i;
}

uint64_t hts_pack_avx2_symbols(const uint8_t *data, uint64_t len, uint8_t *out,
                                  const uint8_t *symbols, int n, int val_per_byte) {
    uint64_t i;
    for (i = 0; i + 32 <= len; i += 32) {
        __m256i bytes = _mm256_loadu_si256((const __m256i *)(data + i));
        __m256i codes = _mm256_setzero_si256();
        for (int s = 1; s < n; ++s) {
            __m256i equal = _mm256_cmpeq_epi8(bytes, _mm256_set1_epi8(symbols[s]));
            codes = _mm256_or_si256(codes, _mm256_and_si256(equal, _mm256_set1_epi8(s)));
        }
        if (val_per_byte == 2) {
            __m256i pairs = _mm256_and_si256(_mm256_or_si256(codes, _mm256_srli_epi16(codes, 4)), _mm256_set1_epi16(255));
            __m256i packed = _mm256_packus_epi16(pairs, _mm256_setzero_si256());
            packed = _mm256_permute4x64_epi64(packed, 0xd8);
            _mm_storeu_si128((__m128i *)(out + i / 2), _mm256_castsi256_si128(packed));
        } else if (val_per_byte == 4) {
            __m256i pairs = _mm256_maddubs_epi16(codes, _mm256_set1_epi16(0x0401));
            __m256i groups = _mm256_madd_epi16(pairs, _mm256_set1_epi32(0x00100001));
            __m256i packed = _mm256_packus_epi32(groups, _mm256_setzero_si256());
            packed = _mm256_packus_epi16(packed, _mm256_setzero_si256());
            __m128i low = _mm256_castsi256_si128(packed);
            __m128i high = _mm256_extracti128_si256(packed, 1);
            _mm_storel_epi64((__m128i *)(out + i / 4), _mm_unpacklo_epi32(low, high));
        } else {
            uint32_t bits = _mm256_movemask_epi8(_mm256_cmpgt_epi8(codes, _mm256_setzero_si256()));
            memcpy(out + i / 8, &bits, sizeof(bits));
        }
    }
    return i;
}

