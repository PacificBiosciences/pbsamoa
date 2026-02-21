Load shared CRAM write helpers and prepare fixtures:

  $ . "${TESTDIR}"/cram_write_helpers.sh
  $ prepare_source_fixtures

fqz slice-size matrix (WCOR-01) validates records-per-slice 1, 10, and 1000:

  $ run_fqz_slice_case 1000
  method-presence-check:* (glob)
  fqz-slice-case:* (glob)
