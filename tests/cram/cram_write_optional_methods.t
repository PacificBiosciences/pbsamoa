Load shared CRAM write helpers and prepare fixtures:

  $ . "${TESTDIR}"/cram_write_helpers.sh
  $ prepare_source_fixtures

Optional write-method interop checks are gated by tool support:

  $ run_optional_write_method_case write_bzip2_optional bzip2 "bzip2"
  method-presence-check:* (glob)
  write-method-case:* (glob)

  $ run_optional_write_method_case write_lzma_optional lzma "lzma"
  method-presence-check:* (glob)
  write-method-case:* (glob)
