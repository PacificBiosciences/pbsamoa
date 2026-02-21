if cram_script.found()
    cram_tests = [
        'bai_build',
        'bai_query',
        'bam_read',
        'bam_write',
        'chunk',
        'edge_cases',
        'format_convert',
        'pipeline',
        'roundtrip',
        'sam_read',
        'sam_write',
        'zmi',
    ]

    cram_env = environment()
    cram_env.prepend('PATH', meson.project_build_root() / 'src' / 'tools')
    cram_env.set('PBSAMOA', pbsamoa.full_path())

    foreach t : cram_tests
        xunit_file_name = 'pbsamoa-cram-' + t.underscorify() + '-test.xml'
        xunit_path = meson.global_build_root() / xunit_file_name

        test(
            'cram-' + t,
            cram_script,
            args: [
                # '--keep-tmpdir',
                '--shell=sh',
                '--verbose',
                '--xunit-file=' + xunit_path,
                files(t + '.t'),
            ],
            depends: [pbsamoa],
            env: cram_env,
            timeout: 36000,
        )
    endforeach
endif
