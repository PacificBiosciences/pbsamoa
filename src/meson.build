################
# version info #
################

gen_headers = [
    vcs_tag(
        command: ['git', 'describe', '--always', '--dirty=*'],
        input: 'LibraryGitHash.hpp.in',
        output: 'LibraryGitHash.hpp',
        replace_string: '@PBSAMOA_GIT_SHA1@',
    ),
    configure_file(
        configuration: {'PBSAMOA_VERSION': meson.project_version()},
        input: 'LibraryVersion.hpp.in',
        output: 'LibraryVersion.hpp',
    ),
]

###########
# library #
###########

lib_install = (not meson.is_subproject()) or (
    get_option('default_library') != 'static'
)

lib = library(
    'pbsamoa',
    files(
        'BaiIndex.cpp',
        'BamRawReader.cpp',
        'BamRecord.cpp',
        'BamRecordReader.cpp',
        'BamWriter.cpp',
        'BamZmwReader.cpp',
        'Bgzf.cpp',
        'CigarClipping.cpp',
        'CigarOp.cpp',
        'LibraryInfo.cpp',
        'RawRecord.cpp',
        'SamHeader.cpp',
        'SamReader.cpp',
        'SamWriter.cpp',
        'Sequence.cpp',
        'TagClipping.cpp',
        'Tags.cpp',
        'ZmiBamWriter.cpp',
        'ZmiWriter.cpp',
        'ZmwIndex.cpp',
        'ZmwWhitelist.cpp',
    ) + gen_headers,
    cpp_args: cpp_flags,
    dependencies: deps,
    include_directories: [inc_dir, thirdparty_dir],
    install: lib_install,
    soversion: meson.project_version(),
    version: meson.project_version(),
)

###############
# executables #
###############

subdir('tools')
