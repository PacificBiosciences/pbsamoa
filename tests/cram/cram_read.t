Define reusable helpers for CRAM method parity checks:

  $ if samtools cram-size --help >/dev/null 2>&1; then HAS_CRAM_SIZE=1; else HAS_CRAM_SIZE=0; fi
  $ normalize_no_ref_sam() {
  >   input_sam="$1"
  >   output_sam="$2"
  >   grep -v "^@PG" "${input_sam}" | awk 'BEGIN{OFS="\\t"} /^@/{print;next} {$7="*";$8=0;$9=0; print}' > "${output_sam}"
  > }
  $ assert_method_presence() {
  >   cram_file="$1"
  >   case_name="$2"
  >   method_name="$3"
  >   method_pattern="$4"
  >   if [ "${HAS_CRAM_SIZE}" -eq 1 ]; then
  >     method_report="${case_name}_cram_size.txt"
  >     samtools cram-size -v "${cram_file}" > "${method_report}"
  >     if grep -E "${method_pattern}" "${method_report}" >/dev/null; then
  >       echo "method-presence-check: verified ${method_name} in ${cram_file}"
  >     else
  >       echo "expected method ${method_name} not present in ${cram_file}" >&2
  >       cat "${method_report}" >&2
  >       return 1
  >     fi
  >   else
  >     echo "method-presence-check: TODO cram-size unavailable, skipped ${method_name} assertion for ${cram_file}"
  >   fi
  > }
  $ run_method_case() {
  >   case_name="$1"
  >   output_fmt="$2"
  >   method_name="$3"
  >   method_pattern="$4"
  >   cram_file="${case_name}.cram"
  >   pbsamoa_sam="${case_name}_pbsamoa.sam"
  >   samtools_sam="${case_name}_samtools.sam"
  >   pbsamoa_norm="${case_name}_pbsamoa_norm.sam"
  >   samtools_norm="${case_name}_samtools_norm.sam"
  >   samtools view --output-fmt "${output_fmt}" -o "${cram_file}" "${TESTDIR}"/../data/many_records.bam
  >   assert_method_presence "${cram_file}" "${case_name}" "${method_name}" "${method_pattern}"
  >   "${PBSAMOA}" dump "${cram_file}" > "${pbsamoa_sam}"
  >   samtools view -h --no-PG "${cram_file}" > "${samtools_sam}"
  >   normalize_no_ref_sam "${pbsamoa_sam}" "${pbsamoa_norm}"
  >   normalize_no_ref_sam "${samtools_sam}" "${samtools_norm}"
  >   diff "${pbsamoa_norm}" "${samtools_norm}"
  > }

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

Read reference-based CRAM using explicit reference FASTA:

  $ printf ">ref\nACGTACGTACGTACGTACGT\n" > ref.fa
  $ printf "@HD\tVN:1.6\tSO:unknown\n@SQ\tSN:ref\tLN:20\nread1\t0\tref\t1\t60\t10M\t*\t0\t0\tACGTTCGTAC\tFFFFFFFFFF\n" > ref.sam
  $ samtools view -bS ref.sam > ref.bam
  $ samtools view -T ref.fa --output-fmt cram,use_rans=0,use_bzip2=0,use_lzma=0,use_fqz=0,use_tok=0,use_arith=0 -o ref.cram ref.bam
  $ "${PBSAMOA}" dump --reference ref.fa ref.cram > pbsamoa_ref.sam
  $ grep -v "^@PG" pbsamoa_ref.sam > pbsamoa_ref_nopg.sam
  $ sed -E 's/\tMD:Z:[^\t]*//g; s/\tNM:i:[^\t]*//g' pbsamoa_ref_nopg.sam > pbsamoa_ref_norm.sam
  $ samtools view -h --no-PG -T ref.fa ref.cram > samtools_ref.sam
  $ grep -v "^@PG" samtools_ref.sam > samtools_ref_nopg.sam
  $ sed -E 's/\tMD:Z:[^\t]*//g; s/\tNM:i:[^\t]*//g' samtools_ref_nopg.sam > samtools_ref_norm.sam
  $ diff pbsamoa_ref_norm.sam samtools_ref_norm.sam

Read reference-based CRAM with wrong FASTA fails with MD5 mismatch:

  $ printf ">ref\nTTTTTTTTTTTTTTTTTTTT\n" > wrong.fa
  $ { "${PBSAMOA}" dump --reference wrong.fa ref.cram > pbsamoa_wrong_ref.sam 2> pbsamoa_wrong_ref.err; dump_status=$?; [ "${dump_status}" -ne 0 ]; } 2> /dev/null
  $ grep "reference MD5 mismatch" pbsamoa_wrong_ref.err
  *reference MD5 mismatch* (glob)

Read reference-free CRAM without providing a FASTA:

  $ samtools view --output-fmt cram,no_ref,use_rans=0,use_bzip2=0,use_lzma=0,use_fqz=0,use_tok=0,use_arith=0 -o ref_no_ref.cram ref.bam
  $ "${PBSAMOA}" dump ref_no_ref.cram > pbsamoa_ref_no_ref.sam
  $ grep -v "^@PG" pbsamoa_ref_no_ref.sam > pbsamoa_ref_no_ref_nopg.sam
  $ awk 'BEGIN{OFS="\\t"} /^@/{print;next} {$7="*";$8=0;$9=0; print}' pbsamoa_ref_no_ref_nopg.sam > pbsamoa_ref_no_ref_norm.sam
  $ samtools view -h --no-PG ref_no_ref.cram > samtools_ref_no_ref.sam
  $ grep -v "^@PG" samtools_ref_no_ref.sam > samtools_ref_no_ref_nopg.sam
  $ awk 'BEGIN{OFS="\\t"} /^@/{print;next} {$7="*";$8=0;$9=0; print}' samtools_ref_no_ref_nopg.sam > samtools_ref_no_ref_norm.sam
  $ diff pbsamoa_ref_no_ref_norm.sam samtools_ref_no_ref_norm.sam

pbsamoa CRAM round-trip preserves BAM dump output:

  $ "${PBSAMOA}" convert "${TESTDIR}"/../data/diverse.bam rt.cram
  $ "${PBSAMOA}" dump rt.cram > from_cram.sam
  $ "${PBSAMOA}" dump "${TESTDIR}"/../data/diverse.bam > from_bam.sam
  $ diff from_cram.sam from_bam.sam

Header SQ lines preserved:

  $ "${PBSAMOA}" dump rt.cram | grep "^@SQ" | wc -l | sed 's/ //g'
  2

Header RG lines preserved:

  $ "${PBSAMOA}" dump rt.cram | grep "^@RG" | wc -l | sed 's/ //g'
  1

Record count preserved:

  $ "${PBSAMOA}" dump rt.cram | grep -cv "^@"
  13

CRAI-backed region queries in pbsamoa dump:

  $ cat > query_region.sam <<'EOF'
  > @HD	VN:1.6	SO:unknown
  > @SQ	SN:chr1	LN:1000
  > region_mapped_1	0	chr1	11	60	5M	*	0	0	ACGTA	IIIII
  > region_mapped_2	0	chr1	19	60	5M	*	0	0	TTGCA	IIIII
  > region_mapped_3	0	chr1	41	60	5M	*	0	0	GGGGG	IIIII
  > region_unmapped_1	4	*	0	0	*	*	0	0	CCCCC	IIIII
  > EOF
  $ samtools view -bS query_region.sam > query_region.bam
  $ "${PBSAMOA}" convert --write-crai query_region.bam query_region.cram
  $ "${PBSAMOA}" dump --region chr1:11-20 query_region.cram | grep -v "^@" | cut -f1
  region_mapped_1
  region_mapped_2
  $ "${PBSAMOA}" dump --index query_region.cram.crai --region '*' query_region.cram | grep -v "^@" | cut -f1
  region_unmapped_1
  $ "${PBSAMOA}" dump --region chr1:1-5 query_region.cram | awk 'BEGIN{n=0} !/^@/{++n} END{print n}'
  0
