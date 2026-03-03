pbsamoa CRAM round-trip preserves BAM dump output:

  $ "${PBSAMOA}" convert "${TESTDIR}"/../data/diverse.bam rt.cram
  $ "${PBSAMOA}" dump rt.cram > from_cram.sam
  $ "${PBSAMOA}" dump "${TESTDIR}"/../data/diverse.bam > from_bam.sam
  $ diff from_cram.sam from_bam.sam

Header SQ lines preserved:

  $ "${PBSAMOA}" dump rt.cram | grep "^@SQ" | wc -l | sed 's/ //g'
  2

Header RG lines preserved:

  $ "${PBSAMOA}" dump rt.cram | grep "^@RG" | wc -l | sed 's/ //g'
  1

Record count preserved:

  $ "${PBSAMOA}" dump rt.cram | grep -cv "^@"
  13
