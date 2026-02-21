# find GoogleTest
gtest_dep = dependency(
    'gtest',
    main: true,
    fallback: ['gtest', 'gtest_main_dep'],
)

# unit test runner
test_exe = executable(
    'test_exe',
    files(
        'TestBaiIndex.cpp',
        'TestBamRawReader.cpp',
        'TestBamRecord.cpp',
        'TestBamRecordClipping.cpp',
        'TestBamRecordReader.cpp',
        'TestBamWriter.cpp',
        'TestBamZmwReader.cpp',
        'TestBgzf.cpp',
        'TestBgzfPipeline.cpp',
        'TestBgzfReader.cpp',
        'TestBgzfWriter.cpp',
        'TestCigarClipping.cpp',
        'TestCigarOp.cpp',
        'TestMetrics.cpp',
        'TestPulseBitset.cpp',
        'TestRawRecord.cpp',
        'TestRawRecordBatch.cpp',
        'TestReaderTypes.cpp',
        'TestSamHeader.cpp',
        'TestSamReader.cpp',
        'TestSamWriter.cpp',
        'TestSequence.cpp',
        'TestSmoke.cpp',
        'TestTagClipping.cpp',
        'TestTags.cpp',
        'TestVirtualOffset.cpp',
        'TestZmi.cpp',
        'TestZmwUtils.cpp',
        'TestZmwWhitelist.cpp',
    ) + configure_file(
        configuration: {
            'TESTS_DATA_DIR': meson.current_source_dir() / '..' / 'data',
        },
        input: 'TestData.hpp.in',
        output: 'TestData.hpp',
    ),
    cpp_args: cpp_flags,
    dependencies: [gtest_dep, deps],
    include_directories: inc_dir,
    link_with: lib,
)

test(
    'unit',
    test_exe,
    protocol: 'gtest',
    timeout: 3600,
)
