Full round-trip: BAM -> SAM -> BAM -> SAM should produce identical non-PG output:

  $ "${PBSAMOA}" dump "${TESTDIR}"/../data/diverse.bam > step1.sam
  $ samtools view -bS step1.sam | samtools sort -o step2.bam
  $ "${PBSAMOA}" dump step2.bam > step3.sam
  $ grep -v "^@PG" step1.sam > step1_nopg.txt
  $ grep -v "^@PG" step3.sam > step3_nopg.txt
  $ diff step1_nopg.txt step3_nopg.txt

Round-trip: SAM -> CRAM -> SAM -> CRAM remains stable with pbsamoa + samtools:

  $ "${PBSAMOA}" convert "${TESTDIR}"/../data/diverse.sam sam_to_cram_step1.cram
  $ samtools view -h --no-PG sam_to_cram_step1.cram > sam_to_cram_step2.sam
  $ "${PBSAMOA}" convert sam_to_cram_step2.sam sam_to_cram_step3.cram
  $ samtools view -h --no-PG sam_to_cram_step3.cram > sam_to_cram_step4.sam
  $ diff sam_to_cram_step2.sam sam_to_cram_step4.sam

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
