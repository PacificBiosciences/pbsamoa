project('isa-l', 'c', version: '2.31.1', license: 'BSD-3-Clause')

add_languages('nasm', native: false)
c = meson.get_compiler('c')
isal_c_args = c.get_supported_arguments(
    '-Wno-sign-compare',
    '-Wno-old-style-declaration',
    '-Wno-implicit-fallthrough',
    '-Wno-unused-parameter',
    '-Wno-pedantic',
)

# Match upstream's igzip build, with portable CRCs and dispatch capped at AVX2.
isal_lib = static_library(
    'pbsamoa_igzip',
    files(
        'igzip/igzip.c',
        'igzip/hufftables_c.c',
        'igzip/igzip_base.c',
        'igzip/igzip_icf_base.c',
        'igzip/adler32_base.c',
        'igzip/flatten_ll.c',
        'igzip/encode_df.c',
        'igzip/igzip_icf_body.c',
        'igzip/huff_codes.c',
        'igzip/igzip_inflate.c',
        'crc/crc_base.c',
        'crc/crc64_base.c',
        'crc/crc_base_aliases.c',
        'igzip/igzip_body.asm',
        'igzip/igzip_finish.asm',
        'igzip/igzip_icf_body_h1_gr_bt.asm',
        'igzip/igzip_icf_finish.asm',
        'igzip/rfc1951_lookup.asm',
        'igzip/adler32_sse.asm',
        'igzip/adler32_avx2_4.asm',
        'igzip/igzip_multibinary.asm',
        'igzip/igzip_update_histogram_01.asm',
        'igzip/igzip_update_histogram_04.asm',
        'igzip/igzip_decode_block_stateless_01.asm',
        'igzip/igzip_decode_block_stateless_04.asm',
        'igzip/igzip_inflate_multibinary.asm',
        'igzip/encode_df_04.asm',
        'igzip/encode_df_06.asm',
        'igzip/proc_heap.asm',
        'igzip/igzip_deflate_hash.asm',
        'igzip/igzip_gen_icf_map_lh1_06.asm',
        'igzip/igzip_gen_icf_map_lh1_04.asm',
        'igzip/igzip_set_long_icf_fg_04.asm',
        'igzip/igzip_set_long_icf_fg_06.asm',
    ),
    c_args: ['-DAS_FEATURE_LEVEL=4'] + isal_c_args,
    nasm_args: ['-DAS_FEATURE_LEVEL=4'],
    include_directories: include_directories('include', 'igzip', 'crc'),
    override_options: ['b_lto=false', 'unity=off', 'c_std=gnu11'],
    install: false,
)

isal_dep = declare_dependency(
    include_directories: include_directories('include', is_system: true),
    link_with: isal_lib,
)
