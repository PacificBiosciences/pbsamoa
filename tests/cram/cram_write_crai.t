Load shared CRAM write helpers and prepare fixtures:

  $ . "${TESTDIR}"/cram_write_helpers.sh
  $ prepare_source_fixtures

CRAI interop matrix (CRAI-04) validates explicit index consumption via samtools `view -X`:

  $ run_write_crai_case write_crai_interop "chr1:1-500"
  crai-interop-case:* (glob)
