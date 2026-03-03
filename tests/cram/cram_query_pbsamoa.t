Load shared CRAM query helpers and build pbsamoa fixtures:

  $ . "${TESTDIR}"/cram_query_helpers.sh
  $ build_pbsamoa_fixture many_pbsamoa "${TESTDIR}"/../data/many_records.bam "100"
  $ build_pbsamoa_fixture diverse_pbsamoa "${TESTDIR}"/../data/diverse.bam ""

pbsamoa-generated CRAI consumed by samtools view -X:

  $ run_parity_case pbsamoa_single_slice_mapped many_pbsamoa.cram many_pbsamoa.crai chr1:1-30 1
  query-parity-case:* (glob)

  $ run_parity_case pbsamoa_multi_slice_boundary many_pbsamoa.cram many_pbsamoa.crai chr1:10050-10060 2
  query-parity-case:* (glob)

  $ run_parity_case pbsamoa_empty_mapped many_pbsamoa.cram many_pbsamoa.crai chr1:60000-60010 0
  query-parity-case:* (glob)

  $ run_parity_case pbsamoa_unmapped diverse_pbsamoa.cram diverse_pbsamoa.crai "*" 1
  query-parity-case:* (glob)

  $ run_parity_case pbsamoa_multi_reference diverse_pbsamoa.cram diverse_pbsamoa.crai chr2:1-500 2
  query-parity-case:* (glob)
