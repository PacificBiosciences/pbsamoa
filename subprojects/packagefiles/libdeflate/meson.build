project('libdeflate', 'c',
  version : '1.22',
  default_options : ['c_std=c11'],
)

libdeflate_sources = files(
  'lib/adler32.c',
  'lib/crc32.c',
  'lib/deflate_compress.c',
  'lib/deflate_decompress.c',
  'lib/gzip_compress.c',
  'lib/gzip_decompress.c',
  'lib/utils.c',
  'lib/zlib_compress.c',
  'lib/zlib_decompress.c',
)

host_machine_cpu = host_machine.cpu_family()
if host_machine_cpu == 'x86_64' or host_machine_cpu == 'x86'
  libdeflate_sources += files('lib/x86/cpu_features.c')
elif host_machine_cpu == 'aarch64'
  libdeflate_sources += files('lib/arm/cpu_features.c')
endif

libdeflate_inc = include_directories('.')

libdeflate_lib = static_library('deflate',
  libdeflate_sources,
  include_directories : libdeflate_inc,
)

libdeflate_dep = declare_dependency(
  include_directories : libdeflate_inc,
  link_with : libdeflate_lib,
)
