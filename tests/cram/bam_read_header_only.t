Test header-only BAM produces only header:

  $ "${PBSAMOA}" dump "${TESTDIR}"/../data/header_only.bam | grep -c "^@"
  * (glob)

  $ "${PBSAMOA}" dump "${TESTDIR}"/../data/header_only.bam | grep -cv "^@"
  0
  [1]
