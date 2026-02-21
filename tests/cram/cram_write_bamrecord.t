Load shared CRAM write helpers and prepare fixtures:

  $ . "${TESTDIR}"/cram_write_helpers.sh
  $ prepare_source_fixtures

Opt-in decode path (`--convert-to-bam-record`) remains parity-correct:

  $ run_write_bamrecord_optin_case write_bamrecord_optin
  bamrecord-optin-case:* (glob)
