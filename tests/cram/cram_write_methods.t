Load shared CRAM write helpers and prepare fixtures:

  $ . "${TESTDIR}"/cram_write_helpers.sh
  $ prepare_source_fixtures

Write-method interop matrix (TEST-02) validates required pbsamoa block methods:

  $ run_write_method_case write_raw raw "raw"
  method-presence-check:* (glob)

  $ run_write_method_case write_gzip gzip "gzip"
  method-presence-check:* (glob)

  $ run_write_method_case write_rans4x8 rans4x8 "r4x8"
  method-presence-check:* (glob)

  $ run_write_method_case write_rans4x16 rans4x16 "rNx16|r4x16"
  method-presence-check:* (glob)

  $ run_write_method_case write_arith arith "arith"
  method-presence-check:* (glob)

  $ run_write_method_case write_fqzcomp fqzcomp "fqzcomp"
  method-presence-check:* (glob)

  $ run_write_method_case write_tok tok "tok3|tok"
  method-presence-check:* (glob)
