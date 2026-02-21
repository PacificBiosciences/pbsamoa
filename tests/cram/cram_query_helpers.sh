#!/usr/bin/env sh

# Shared helpers for CRAM query tests.
# Sourced by each cram_query_*.t file.

normalize_query_records() {
  input_sam="$1"
  output_sam="$2"
  awk '!/^@/{print}' "${input_sam}" | LC_ALL=C sort > "${output_sam}"
}

build_pbsamoa_fixture() {
  fixture_name="$1"
  input_bam="$2"
  records_per_slice="$3"
  cram_file="${fixture_name}.cram"
  crai_file="${fixture_name}.crai"
  if [ -n "${records_per_slice}" ]; then
    "${PBSAMOA}" convert --write-crai --records-per-slice "${records_per_slice}" "${input_bam}" "${cram_file}"
  else
    "${PBSAMOA}" convert --write-crai "${input_bam}" "${cram_file}"
  fi
  [ -f "${cram_file}.crai" ]
  cp "${cram_file}.crai" "${crai_file}"
  rm "${cram_file}.crai"
}

build_samtools_fixture() {
  fixture_name="$1"
  input_bam="$2"
  cram_file="${fixture_name}.cram"
  crai_file="${fixture_name}.crai"
  samtools view --output-fmt cram,no_ref -o "${cram_file}" "${input_bam}"
  samtools index "${cram_file}"
  [ -f "${cram_file}.crai" ]
  cp "${cram_file}.crai" "${crai_file}"
  rm "${cram_file}.crai"
}

run_parity_case() {
  case_name="$1"
  cram_file="$2"
  crai_file="$3"
  region="$4"
  expected_count="$5"
  pbsamoa_raw="${case_name}_pbsamoa_raw.sam"
  pbsamoa_norm="${case_name}_pbsamoa_norm.sam"
  samtools_raw="${case_name}_samtools_raw.sam"
  samtools_norm="${case_name}_samtools_norm.sam"
  "${PBSAMOA}" dump --index "${crai_file}" --region "${region}" "${cram_file}" > "${pbsamoa_raw}"
  samtools view -X "${cram_file}" "${crai_file}" "${region}" > "${samtools_raw}"
  normalize_query_records "${pbsamoa_raw}" "${pbsamoa_norm}"
  LC_ALL=C sort "${samtools_raw}" > "${samtools_norm}"
  diff "${pbsamoa_norm}" "${samtools_norm}"
  pbsamoa_count="$(awk 'END{print NR+0}' "${pbsamoa_norm}")"
  samtools_count="$(awk 'END{print NR+0}' "${samtools_norm}")"
  [ "${pbsamoa_count}" -eq "${expected_count}" ]
  [ "${samtools_count}" -eq "${expected_count}" ]
  echo "query-parity-case: ${case_name} region ${region} count ${expected_count}"
}
