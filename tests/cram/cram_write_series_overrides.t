Load shared CRAM write helpers and prepare fixtures:

  $ . "${TESTDIR}"/cram_write_helpers.sh
  $ prepare_source_fixtures

Per-series compression overrides apply mixed codecs within a single CRAM:

  $ run_series_override_case write_series_overrides
  method-presence-check:* (glob)
  method-presence-check:* (glob)
  series-override-case:* (glob)
