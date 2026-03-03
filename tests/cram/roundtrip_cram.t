Round-trip: SAM -> CRAM -> SAM -> CRAM remains stable with pbsamoa + samtools:

  $ "${PBSAMOA}" convert "${TESTDIR}"/../data/diverse.sam sam_to_cram_step1.cram
  $ samtools view -h --no-PG sam_to_cram_step1.cram > sam_to_cram_step2.sam
  $ "${PBSAMOA}" convert sam_to_cram_step2.sam sam_to_cram_step3.cram
  $ samtools view -h --no-PG sam_to_cram_step3.cram > sam_to_cram_step4.sam
  $ diff sam_to_cram_step2.sam sam_to_cram_step4.sam
