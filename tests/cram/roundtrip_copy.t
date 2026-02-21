Copy via pbsamoa zmi-build preserves all records:

  $ "${PBSAMOA}" zmi-build "${TESTDIR}"/../data/diverse.bam copy.bam 2>/dev/null
  $ "${PBSAMOA}" dump copy.bam > copy.sam
  $ "${PBSAMOA}" dump "${TESTDIR}"/../data/diverse.bam > original.sam
  $ diff copy.sam original.sam

Record count is preserved:

  $ "${PBSAMOA}" dump "${TESTDIR}"/../data/diverse.bam | grep -cv "^@"
  13

Header SQ lines preserved:

  $ "${PBSAMOA}" dump "${TESTDIR}"/../data/diverse.bam | grep "^@SQ" | wc -l | sed 's/ //g'
  2

Header RG lines preserved:

  $ "${PBSAMOA}" dump "${TESTDIR}"/../data/diverse.bam | grep "^@RG" | wc -l | sed 's/ //g'
  1
