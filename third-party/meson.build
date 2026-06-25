thirdparty_dir = include_directories('.')
thirdparty_htscodecs_dir = include_directories('htscodecs')
thirdparty_c_args = c.get_supported_arguments(
    # Imported upstream C sources are warning-heavy under strict project flags.
    '-Wno-sign-compare',
    '-Wno-unused-parameter',
    # GCC 16: dead 'nlines' loop counter in htscodecs tokenise_name3.c
    '-Wno-unused-but-set-variable',
)

htscodecs_sources = files(
    'htscodecs/arith_dynamic.c',
    'htscodecs/fqzcomp_qual.c',
    'htscodecs/pack.c',
    'htscodecs/rANS_static.c',
    'htscodecs/rANS_static4x16pr.c',
    'htscodecs/rANS_static32x16pr.c',
    'htscodecs/rANS_static32x16pr_avx2.c',
    'htscodecs/rANS_static32x16pr_avx512.c',
    'htscodecs/rANS_static32x16pr_neon.c',
    'htscodecs/rANS_static32x16pr_sse4.c',
    'htscodecs/rle.c',
    'htscodecs/tokenise_name3.c',
    'htscodecs/utils.c',
)

htscodecs_lib = static_library(
    'pbsamoa_htscodecs',
    htscodecs_sources,
    c_args: thirdparty_c_args,
    include_directories: [thirdparty_dir, thirdparty_htscodecs_dir],
    # These C sources intentionally reuse macro names between files and are
    # not unity-safe.
    override_options: ['unity=off'],
    install: false,
)

if meson.is_subproject()
    subdir_done()
endif
