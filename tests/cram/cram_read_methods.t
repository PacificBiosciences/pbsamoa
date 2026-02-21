Load shared CRAM read helpers:

  $ . "${TESTDIR}"/cram_read_helpers.sh

Read CRAM with gzip blocks created by samtools:

  $ run_method_case gzip_case "cram,no_ref,use_rans=0,use_bzip2=0,use_lzma=0,use_fqz=0,use_tok=0,use_arith=0" "gzip" "gzip"
  method-presence-check:* (glob)

Read CRAM with bzip2 blocks created by samtools:

  $ run_method_case bzip2_case "cram,no_ref,use_rans=0,use_bzip2=1,use_lzma=0,use_fqz=0,use_tok=0,use_arith=0" "bzip2" "bzip2"
  method-presence-check:* (glob)

Read CRAM with lzma blocks created by samtools:

  $ run_method_case lzma_case "cram,no_ref,use_rans=0,use_bzip2=0,use_lzma=1,use_fqz=0,use_tok=0,use_arith=0" "lzma" "lzma"
  method-presence-check:* (glob)

Read CRAM with rANS 4x8 blocks created by samtools:

  $ run_method_case rans4x8_case "cram,version=3.0,no_ref,use_rans=1,use_bzip2=0,use_lzma=0,use_fqz=0,use_tok=0,use_arith=0" "rANS 4x8" "r4x8"
  method-presence-check:* (glob)

Read CRAM with rANS 4x16 blocks created by samtools:

  $ run_method_case rans4x16_case "cram,version=3.1,no_ref,use_rans=1,use_bzip2=0,use_lzma=0,use_fqz=0,use_tok=0,use_arith=0" "rANS 4x16" "rNx16|r4x16"
  method-presence-check:* (glob)

Read CRAM with adaptive arith blocks created by samtools:

  $ run_method_case arith_case "cram,version=3.1,no_ref,use_rans=0,use_bzip2=0,use_lzma=0,use_fqz=0,use_tok=0,use_arith=1" "adaptive arith" "arith"
  method-presence-check:* (glob)

Read CRAM with fqzcomp blocks created by samtools:

  $ run_method_case fqz_case "cram,version=3.1,no_ref,use_rans=0,use_bzip2=0,use_lzma=0,use_fqz=1,use_tok=0,use_arith=0" "fqzcomp" "fqzcomp"
  method-presence-check:* (glob)

Read CRAM with name tokeniser blocks created by samtools:

  $ run_method_case tok_case "cram,version=3.1,no_ref,use_rans=0,use_bzip2=0,use_lzma=0,use_fqz=0,use_tok=1,use_arith=0" "name tokeniser" "tok3|tok"
  method-presence-check:* (glob)
