#!/bin/bash
set -euo pipefail

SCRIPT_DIR=$(dirname "$(realpath "$0")")

main() {
  local -r input_data_dir="/groups/umcg-gdio/tmp02/test/FBD12345/fastq/"
  local -r input_adaptive_sampling_csv="./adaptive_sampling_FBD12345_f3dbc49_e7687v7.csv"

  local -r output_dir="${SCRIPT_DIR}/output"
  local -r output_intermediates_dir="${output_dir}/intermediates"
  
  if [[ -d "${output_dir}" ]]; then
   echo -e "error: output directory '${output_dir}' already exists"
   exit 1
  fi
  if [[ -d "${output_intermediates_dir}" ]]; then
   echo -e "error: output intermediates directory '${output_intermediates_dir}' already exists"
   exit 1
  fi

  mkdir -p "${output_dir}"
  mkdir -p "${output_intermediates_dir}"

  wget --directory-prefix "${output_intermediates_dir}" "https://downloads.molgeniscloud.org/downloads/vip/_dev/images/utils/seqtk-1.4.sif"

  #stap 1: merge and decompress the passed fastq files
  cat "${input_data_dir}"/*.fastq.gz > "${output_intermediates_dir}/merged.fastq.gz"
  gunzip "${output_intermediates_dir}/merged.fastq.gz"

  #stap 2: extract read IDs of interest based on decision in csv file
  grep stop_receiving "${input_adaptive_sampling_csv}" | cut -d , -f 5 > "${output_intermediates_dir}/read_IDs.txt"

  #stap 3: extract reads and save to new file
  APPTAINER_BIND=/groups apptainer exec "${output_intermediates_dir}/seqtk-1.4.sif" seqtk subseq "${output_intermediates_dir}/merged.fastq" "${output_intermediates_dir}/read_IDs.txt" > "${output_intermediates_dir}/out_seqtk.fastq"

  gzip "${output_intermediates_dir}/out_seqtk.fastq"
}

main "${@}"
