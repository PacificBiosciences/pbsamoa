Load shared CRAM query helpers and build samtools fixtures:

  $ . "${TESTDIR}"/cram_query_helpers.sh
  $ build_samtools_fixture many_samtools "${TESTDIR}"/../data/many_records.bam
  $ build_samtools_fixture diverse_samtools "${TESTDIR}"/../data/diverse.bam

samtools-generated CRAI consumed by pbsamoa dump:

  $ run_parity_case samtools_single_slice_mapped many_samtools.cram many_samtools.crai chr1:1-30 1
  query-parity-case:* (glob)

  $ run_parity_case samtools_multi_slice_boundary many_samtools.cram many_samtools.crai chr1:10050-10060 2
  query-parity-case:* (glob)

  $ run_parity_case samtools_empty_mapped many_samtools.cram many_samtools.crai chr1:60000-60010 0
  query-parity-case:* (glob)

  $ run_parity_case samtools_unmapped diverse_samtools.cram diverse_samtools.crai "*" 1
  query-parity-case:* (glob)

  $ run_parity_case samtools_multi_reference diverse_samtools.cram diverse_samtools.crai chr2:1-500 2
  query-parity-case:* (glob)
