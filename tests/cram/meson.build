if cram_script.found()
    cram_tests = [
        'bai_build_create',
        'bai_build_parity',
        'bai_query',
        'bam_read_dump',
        'bam_read_header_only',
        'bam_write_parity',
        'bam_write_zmi',
        'chunk_basic',
        'chunk_five_way',
        'chunk_three_way',
        'chunk_two_way',
        'dump_header_flags',
        'cram_query_pbsamoa',
        'cram_query_samtools',
        'cram_read_crai_query',
        'cram_read_methods',
        'cram_read_reference',
        'cram_read_roundtrip',
        'cram_write_bamrecord',
        'cram_write_crai',
        'cram_write_fqz_slices_1',
        'cram_write_fqz_slices_10',
        'cram_write_fqz_slices_1000',
        'cram_write_methods',
        'cram_write_optional_methods',
        'cram_write_series_overrides',
        'edge_cases_diverse',
        'edge_cases_empty',
        'format_convert_interop',
        'format_convert_roundtrip',
        'pipeline_copy',
        'pipeline_dump',
        'roundtrip_bam',
        'roundtrip_copy',
        'roundtrip_cram',
        'sam_read',
        'sam_write_parity',
        'sam_write_roundtrip',
        'zmi_build',
        'zmi_chunk',
        'zmi_index',
        'zmi_query',
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
