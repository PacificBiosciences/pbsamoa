project(
    'PBSAMOA',
    ['cpp', 'c'],
    default_options: [
        'b_ndebug=if-release',
        'buildtype=release',
        'c_std=c23',
        'cpp_std=c++23',
        'warning_level=3',
    ],
    license: 'BSD-3',
    meson_version: '>= 1.1.0',
    version: '0.1.0',
)

############
# CXXFLAGS #
############

cpp = meson.get_compiler('cpp')
c = meson.get_compiler('c')
cpp_flags = cpp.get_supported_arguments(
    # stupid warning https://gcc.gnu.org/bugzilla/show_bug.cgi?id=96868#c3
    '-Wno-missing-field-initializers',
    # currently broken upstream https://gcc.gnu.org/bugzilla/show_bug.cgi?id=61596
    '-Wno-unused-local-typedefs',
    '-Wdouble-promotion',
    '-Wduplicated-cond',
    '-Wformat=1',
    '-Wlogical-op',
    '-Wold-style-cast',
    '-Wrestrict',
    '-Wshadow',
    '-Wsuggest-override',
    '-Wunused-local-typedef',
    '-Wunused-parameter',
    '-Wuseless-cast',
    # GCC 14 false positive: std::expected<T,E> move with T containing std::string
    # triggers spurious -Wmaybe-uninitialized from standard library internals
    '-Wno-maybe-uninitialized',
)

################
# dependencies #
################

# libdeflate
libdeflate_dep = dependency(
    'libdeflate',
    fallback: ['libdeflate', 'libdeflate_dep'],
    include_type: 'system',
)

deps = [libdeflate_dep, dependency('threads')]

# Native CRAM implementation intentionally does not depend on htslib.
bzip2_dep = dependency('bzip2', required: false, include_type: 'system')
if bzip2_dep.found()
    deps += bzip2_dep
    add_project_arguments('-DPBSAMOA_HAVE_BZIP2=1', language: ['cpp'])
endif

liblzma_dep = dependency('liblzma', required: false, include_type: 'system')
if liblzma_dep.found()
    deps += liblzma_dep
    add_project_arguments('-DPBSAMOA_HAVE_LZMA=1', language: ['cpp'])
endif

libm_dep = c.find_library('m', required: false)
if libm_dep.found()
    deps += libm_dep
endif

###########
# headers #
###########

subdir('include')

###############
# third-party #
###############

subdir('third-party')

##########################
# libraries & executable #
##########################

subdir('src')

#########
# tests #
#########

if (not meson.is_subproject()) and get_option('tests')
    pbtestutils_dep = dependency(
        'pbtestutils',
        fallback: ['pbtestutils', 'pbtestutils_dep'],
    )
    cram_script = subproject('pbtestutils').get_variable('cram_script')
    subdir('tests')
endif

###################
# dependency info #
###################

if not meson.is_subproject()
    import('pkgconfig').generate(
        lib,
        description: 'PacBio SAM/BAM/BAI C++23 library',
        filebase: 'pbsamoa',
        name: 'pbsamoa',
        version: meson.project_version(),
    )
endif

pbsamoa_dep = declare_dependency(
    dependencies: deps,
    include_directories: [inc_dir, thirdparty_dir],
    link_with: lib,
    version: meson.project_version(),
)

meson.override_dependency('pbsamoa', pbsamoa_dep)
