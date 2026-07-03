Error handling: invalid arguments and missing input must exit non-zero. Errors are
reported through the pbcopper CLIv2 logger, i.e. a timestamped FATAL line ending in
"<tool> ERROR: <message>"; the leading (glob) star absorbs the log prefix.

Invalid --order value:

  $ "${PBSAMOA}" sort --order bogus "${TESTDIR}"/../data/unsorted.bam out.bam
  *pbsamoa sort ERROR: invalid --order: bogus (glob)
  [1]

--order tag without --tag:

  $ "${PBSAMOA}" sort --order tag "${TESTDIR}"/../data/unsorted.bam out.bam
  *pbsamoa sort ERROR: --order tag requires --tag XX (glob)
  [1]

Missing input file:

  $ "${PBSAMOA}" sort "${TESTDIR}"/../data/no_such_file.bam out.bam
  *pbsamoa sort ERROR: SortBam: input file does not exist: */no_such_file.bam (glob)
  [1]

Missing output operand prints usage and exits non-zero:

  $ "${PBSAMOA}" sort "${TESTDIR}"/../data/unsorted.bam > /dev/null 2>&1
  [1]
