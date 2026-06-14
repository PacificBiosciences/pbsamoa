Error handling: invalid arguments and missing input must exit non-zero.

Invalid --order value:

  $ "${PBSAMOA}" sort --order bogus "${TESTDIR}"/../data/unsorted.bam out.bam
  Error: invalid --order: bogus
  [1]

--order tag without --tag:

  $ "${PBSAMOA}" sort --order tag "${TESTDIR}"/../data/unsorted.bam out.bam
  Error: --order tag requires --tag XX
  [1]

Missing input file:

  $ "${PBSAMOA}" sort "${TESTDIR}"/../data/no_such_file.bam out.bam
  Error: SortBam: input file does not exist: */no_such_file.bam (glob)
  [1]

Missing output operand prints usage and exits non-zero:

  $ "${PBSAMOA}" sort "${TESTDIR}"/../data/unsorted.bam > /dev/null 2>&1
  [1]
