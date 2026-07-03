`pbsamoa dump` header flags behave like `bam2sam`-style controls:

  $ "${PBSAMOA}" dump --no-header "${TESTDIR}"/../data/spec_example.bam > pbsamoa_no_header.sam
  $ samtools view --no-PG "${TESTDIR}"/../data/spec_example.bam > samtools_no_header.sam
  $ diff pbsamoa_no_header.sam samtools_no_header.sam

  $ "${PBSAMOA}" dump --header-only "${TESTDIR}"/../data/spec_example.bam > pbsamoa_header_only.sam
  $ samtools view -H --no-PG "${TESTDIR}"/../data/spec_example.bam > samtools_header_only.sam
  $ diff pbsamoa_header_only.sam samtools_header_only.sam

  $ "${PBSAMOA}" dump --no-header "${TESTDIR}"/../data/spec_example.sam | awk 'BEGIN{n=0} /^@/{n++} END{print n}'
  0

  $ "${PBSAMOA}" convert "${TESTDIR}"/../data/spec_example.bam dump_header_flags.cram
  $ "${PBSAMOA}" dump --header-only dump_header_flags.cram | awk 'BEGIN{n=0} !/^@/{n++} END{print n}'
  0

  $ "${PBSAMOA}" dump --no-header dump_header_flags.cram | awk 'BEGIN{n=0} /^@/{n++} END{print n}'
  0

  $ "${PBSAMOA}" dump --no-header --header-only "${TESTDIR}"/../data/spec_example.bam > /dev/null
  *--no-header and --header-only are mutually exclusive* (glob)
  [1]
