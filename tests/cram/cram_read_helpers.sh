#!/usr/bin/env sh

# Shared helpers for CRAM read tests.
# Sourced by each cram_read_*.t file.

if samtools cram-size --help >/dev/null 2>&1; then HAS_CRAM_SIZE=1; else HAS_CRAM_SIZE=0; fi

normalize_no_ref_sam() {
  input_sam="$1"
  output_sam="$2"
  grep -v "^@PG" "${input_sam}" | awk 'BEGIN{OFS="\t"} /^@/{print;next} {$7="*";$8=0;$9=0; print}' > "${output_sam}"
}

assert_method_presence() {
  cram_file="$1"
  case_name="$2"
  method_name="$3"
  method_pattern="$4"
  if [ "${HAS_CRAM_SIZE}" -eq 1 ]; then
    method_report="${case_name}_cram_size.txt"
    samtools cram-size -v "${cram_file}" > "${method_report}"
    if grep -E "${method_pattern}" "${method_report}" >/dev/null; then
      echo "method-presence-check: verified ${method_name} in ${cram_file}"
    else
      echo "expected method ${method_name} not present in ${cram_file}" >&2
      cat "${method_report}" >&2
      return 1
    fi
  else
    echo "method-presence-check: TODO cram-size unavailable, skipped ${method_name} assertion for ${cram_file}"
  fi
}

run_method_case() {
  case_name="$1"
  output_fmt="$2"
  method_name="$3"
  method_pattern="$4"
  cram_file="${case_name}.cram"
  pbsamoa_sam="${case_name}_pbsamoa.sam"
  samtools_sam="${case_name}_samtools.sam"
  pbsamoa_norm="${case_name}_pbsamoa_norm.sam"
  samtools_norm="${case_name}_samtools_norm.sam"
  samtools view --output-fmt "${output_fmt}" -o "${cram_file}" "${TESTDIR}"/../data/many_records.bam
  assert_method_presence "${cram_file}" "${case_name}" "${method_name}" "${method_pattern}"
  "${PBSAMOA}" dump "${cram_file}" > "${pbsamoa_sam}"
  samtools view -h --no-PG "${cram_file}" > "${samtools_sam}"
  normalize_no_ref_sam "${pbsamoa_sam}" "${pbsamoa_norm}"
  normalize_no_ref_sam "${samtools_sam}" "${samtools_norm}"
  diff "${pbsamoa_norm}" "${samtools_norm}"
}
