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

Diverse record types all survive round-trip:

  $ "${PBSAMOA}" dump "${TESTDIR}"/../data/diverse.bam | grep "unmapped" | cut -f3
  *

  $ "${PBSAMOA}" dump "${TESTDIR}"/../data/diverse.bam | grep "hard_clip" | cut -f6
  3H7M

  $ "${PBSAMOA}" dump "${TESTDIR}"/../data/diverse.bam | grep "no_qual" | cut -f11
  *

  $ "${PBSAMOA}" dump "${TESTDIR}"/../data/diverse.bam | grep "array_tag" | grep -o "XA:B:i,[0-9,]*"
  XA:B:i,1,2,3,4,5
