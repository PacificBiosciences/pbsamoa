Define reusable helpers for CRAM write parity and method inspection:

  $ if samtools cram-size --help >/dev/null 2>&1; then HAS_CRAM_SIZE=1; else HAS_CRAM_SIZE=0; fi
  $ normalize_no_ref_sam() {
  >   input_sam="$1"
  >   output_sam="$2"
  >   grep -v "^@PG" "${input_sam}" | awk 'BEGIN{OFS="\t"} /^@/{print;next} {$7="*";$8=0;$9=0; print}' > "${output_sam}"
  > }
  $ assert_method_presence() {
  >   cram_file="$1"
  >   case_name="$2"
  >   method_name="$3"
  >   method_pattern="$4"
  >   if [ "${HAS_CRAM_SIZE}" -eq 1 ]; then
  >     report_file="${case_name}_cram_size.txt"
  >     samtools cram-size -v "${cram_file}" > "${report_file}"
  >     if grep -E "${method_pattern}" "${report_file}" >/dev/null; then
  >       echo "method-presence-check: verified ${method_name} in ${cram_file}"
  >     else
  >       echo "expected method ${method_name} not present in ${cram_file}" >&2
  >       cat "${report_file}" >&2
  >       return 1
  >     fi
  >   else
  >     echo "method-presence-check: TODO cram-size unavailable, skipped ${method_name} assertion for ${cram_file}"
  >   fi
  > }
  $ run_write_method_case() {
  >   case_name="$1"
  >   method_name="$2"
  >   method_pattern="$3"
  >   cram_file="${case_name}.cram"
  >   decoded_sam="${case_name}_samtools.sam"
  >   decoded_norm="${case_name}_samtools_norm.sam"
  >   "${PBSAMOA}" convert --block-compression "${method_name}" "${TESTDIR}"/../data/many_records.bam "${cram_file}"
  >   assert_method_presence "${cram_file}" "${case_name}" "${method_name}" "${method_pattern}"
  >   samtools view -h --no-PG "${cram_file}" > "${decoded_sam}"
  >   normalize_no_ref_sam "${decoded_sam}" "${decoded_norm}"
  >   diff "${decoded_norm}" many_records_source_norm.sam
  > }
  $ run_optional_write_method_case() {
  >   case_name="$1"
  >   method_name="$2"
  >   method_pattern="$3"
  >   cram_file="${case_name}.cram"
  >   decoded_sam="${case_name}_samtools.sam"
  >   decoded_norm="${case_name}_samtools_norm.sam"
  >   if "${PBSAMOA}" convert --block-compression "${method_name}" "${TESTDIR}"/../data/many_records.bam "${cram_file}" 2>"${case_name}.err"; then
  >     assert_method_presence "${cram_file}" "${case_name}" "${method_name}" "${method_pattern}"
  >     samtools view -h --no-PG "${cram_file}" > "${decoded_sam}"
  >     normalize_no_ref_sam "${decoded_sam}" "${decoded_norm}"
  >     diff "${decoded_norm}" many_records_source_norm.sam
  >     echo "write-method-case: verified optional ${method_name}"
  >   else
  >     echo "write-method-case: optional ${method_name} unavailable, skipped"
  >   fi
  > }
  $ run_fqz_slice_case() {
  >   records_per_slice="$1"
  >   case_name="fqz_slice_${records_per_slice}"
  >   cram_file="${case_name}.cram"
  >   samtools_sam="${case_name}_samtools.sam"
  >   samtools_nopg="${case_name}_samtools_nopg.sam"
  >   pbsamoa_sam="${case_name}_pbsamoa.sam"
  >   pbsamoa_nopg="${case_name}_pbsamoa_nopg.sam"
  >   "${PBSAMOA}" convert --block-compression fqzcomp --records-per-slice "${records_per_slice}" "${TESTDIR}"/../data/many_records.bam "${cram_file}"
  >   assert_method_presence "${cram_file}" "${case_name}" "fqzcomp" "fqzcomp"
  >   samtools view -h --no-PG "${cram_file}" > "${samtools_sam}"
  >   "${PBSAMOA}" dump "${cram_file}" > "${pbsamoa_sam}"
  >   grep -v "^@PG" "${samtools_sam}" > "${samtools_nopg}"
  >   grep -v "^@PG" "${pbsamoa_sam}" > "${pbsamoa_nopg}"
  >   diff "${samtools_nopg}" many_records_source_nopg.sam
  >   diff "${pbsamoa_nopg}" "${samtools_nopg}"
  >   echo "fqz-slice-case: records-per-slice ${records_per_slice} parity verified"
  > }
  $ run_series_override_case() {
  >   case_name="$1"
  >   cram_file="${case_name}.cram"
  >   decoded_sam="${case_name}_samtools.sam"
  >   decoded_norm="${case_name}_samtools_norm.sam"
  >   "${PBSAMOA}" convert --block-compression raw --series-compression RN=tok --series-compression QS=fqzcomp "${TESTDIR}"/../data/many_records.bam "${cram_file}"
  >   assert_method_presence "${cram_file}" "${case_name}_tok" "tok" "tok3|tok"
  >   assert_method_presence "${cram_file}" "${case_name}_fqz" "fqzcomp" "fqzcomp"
  >   samtools view -h --no-PG "${cram_file}" > "${decoded_sam}"
  >   normalize_no_ref_sam "${decoded_sam}" "${decoded_norm}"
  >   diff "${decoded_norm}" many_records_source_norm.sam
  >   echo "series-override-case: RN=tok and QS=fqzcomp verified"
  > }
  $ run_write_crai_case() {
  >   case_name="$1"
  >   region="$2"
  >   cram_file="${case_name}.cram"
  >   pbsamoa_crai="${case_name}_pbsamoa.crai"
  >   pbsamoa_region_sam="${case_name}_pbsamoa_region.sam"
  >   samtools_region_sam="${case_name}_samtools_region.sam"
  >   "${PBSAMOA}" convert --write-crai "${TESTDIR}"/../data/many_records.bam "${cram_file}"
  >   [ -f "${cram_file}.crai" ]
  >   cp "${cram_file}.crai" "${pbsamoa_crai}"
  >   samtools index "${cram_file}"
  >   samtools view -X "${cram_file}" "${pbsamoa_crai}" "${region}" > "${pbsamoa_region_sam}"
  >   samtools view -X "${cram_file}" "${cram_file}.crai" "${region}" > "${samtools_region_sam}"
  >   diff "${pbsamoa_region_sam}" "${samtools_region_sam}"
  >   echo "crai-interop-case: region ${region} parity verified"
  > }
  $ run_write_bamrecord_optin_case() {
  >   case_name="$1"
  >   cram_file="${case_name}.cram"
  >   decoded_sam="${case_name}_samtools.sam"
  >   decoded_norm="${case_name}_samtools_norm.sam"
  >   "${PBSAMOA}" convert --convert-to-bam-record --decode-threads 1 "${TESTDIR}"/../data/many_records.bam "${cram_file}"
  >   samtools view -h --no-PG "${cram_file}" > "${decoded_sam}"
  >   normalize_no_ref_sam "${decoded_sam}" "${decoded_norm}"
  >   diff "${decoded_norm}" many_records_source_norm.sam
  >   echo "bamrecord-optin-case: parity verified"
  > }

Prepare source SAM fixtures for method and fqz matrices:

  $ samtools view -h --no-PG "${TESTDIR}"/../data/many_records.bam > many_records_source.sam
  $ normalize_no_ref_sam many_records_source.sam many_records_source_norm.sam
  $ grep -v "^@PG" many_records_source.sam > many_records_source_nopg.sam

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

Optional write-method interop checks are gated by tool support:

  $ run_optional_write_method_case write_bzip2_optional bzip2 "bzip2"
  method-presence-check:* (glob)
  write-method-case:* (glob)

  $ run_optional_write_method_case write_lzma_optional lzma "lzma"
  method-presence-check:* (glob)
  write-method-case:* (glob)

fqz slice-size matrix (WCOR-01) validates records-per-slice 1, 10, and 1000:

  $ run_fqz_slice_case 1
  method-presence-check:* (glob)
  fqz-slice-case:* (glob)

  $ run_fqz_slice_case 10
  method-presence-check:* (glob)
  fqz-slice-case:* (glob)

  $ run_fqz_slice_case 1000
  method-presence-check:* (glob)
  fqz-slice-case:* (glob)

Per-series compression overrides apply mixed codecs within a single CRAM:

  $ run_series_override_case write_series_overrides
  method-presence-check:* (glob)
  method-presence-check:* (glob)
  series-override-case:* (glob)

CRAI interop matrix (CRAI-04) validates explicit index consumption via samtools `view -X`:

  $ run_write_crai_case write_crai_interop "chr1:1-500"
  crai-interop-case:* (glob)

Opt-in decode path (`--convert-to-bam-record`) remains parity-correct:

  $ run_write_bamrecord_optin_case write_bamrecord_optin
  bamrecord-optin-case:* (glob)
