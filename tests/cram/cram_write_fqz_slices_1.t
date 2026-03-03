Load shared CRAM write helpers:

  $ . "${TESTDIR}"/cram_write_helpers.sh

fqz slice-size 1 validates single-record-per-slice boundary (WCOR-01):

  $ "${PBSAMOA}" convert --block-compression fqzcomp --compression-threads 8 --records-per-slice 1 "${TESTDIR}"/../data/diverse.bam fqz_slice_1.cram
  $ assert_method_presence fqz_slice_1.cram fqz_slice_1 "fqzcomp" "fqzcomp"
  method-presence-check:* (glob)
  $ "${PBSAMOA}" dump fqz_slice_1.cram > fqz_slice_1_pbsamoa.sam
  $ "${PBSAMOA}" dump "${TESTDIR}"/../data/diverse.bam > fqz_slice_1_source.sam
  $ grep -v "^@PG" fqz_slice_1_pbsamoa.sam > fqz_slice_1_pbsamoa_nopg.sam
  $ grep -v "^@PG" fqz_slice_1_source.sam > fqz_slice_1_source_nopg.sam
  $ diff fqz_slice_1_pbsamoa_nopg.sam fqz_slice_1_source_nopg.sam
