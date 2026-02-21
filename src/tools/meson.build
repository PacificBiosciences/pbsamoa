pbsamoa_tool_sources = files('pbsamoa.cpp')

subdir('bai-build')
subdir('bai-query')
subdir('bench')
subdir('chunk')
subdir('dump')
subdir('zmi-build')
subdir('zmi-query')

pbsamoa = executable(
    'pbsamoa',
    pbsamoa_tool_sources,
    cpp_args: cpp_flags,
    dependencies: deps,
    include_directories: [inc_dir, thirdparty_dir],
    install: not meson.is_subproject(),
    link_with: lib,
)
