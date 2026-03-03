Empty BAM file produces header-only output:

  $ "${PBSAMOA}" dump "${TESTDIR}"/../data/empty_records.bam | grep -c "^@"
  3

  $ "${PBSAMOA}" dump "${TESTDIR}"/../data/empty_records.bam | grep -v "^@" | wc -l | sed 's/ //g'
  0

Header-only BAM:

  $ "${PBSAMOA}" dump "${TESTDIR}"/../data/header_only.bam | grep -c "^@"
  3

  $ "${PBSAMOA}" dump "${TESTDIR}"/../data/header_only.bam | grep -v "^@" | wc -l | sed 's/ //g'
  0
