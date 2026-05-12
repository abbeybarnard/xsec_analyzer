#!/bin/bash

# Usage:
# ./scripts/ReprocessNTuples.sh OUTPUT_DIRECTORY SELECTION_NAMES NTUPLE_LIST_FILE

if [ "$#" -ne 3 ]; then
  echo "Usage: ./scripts/ReprocessNTuples.sh OUTPUT_DIRECTORY SELECTION_NAMES NTUPLE_LIST_FILE"
  exit 1
fi

output_dir=$1
selections=$2
ntuple_list_file=$3

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
process_ntuples_bin="${repo_dir}/bin/ProcessNTuples"

if [ ! -x "${process_ntuples_bin}" ]; then
  echo "Executable \"${process_ntuples_bin}\" not found. Build with make first."
  exit 3
fi

# Ensure output directory exists
mkdir -p "${output_dir}"

# ------------------------------------------
# Write analysis file list
# ------------------------------------------
analysis_list_file="${repo_dir}/configs/files_to_process_analysis.txt"
: > "${analysis_list_file}"

echo "Writing analysis file list to:"
echo "  ${analysis_list_file}"
echo ""

# Count files
total_files=$(grep -v '^\s*#' "$ntuple_list_file" | grep -v '^\s*$' | wc -l)
echo "Total number of files = ${total_files}"

counter=0

# ------------------------------------------
# MAIN LOOP
# ------------------------------------------
while read -r file_path file_type beam_mode; do

  # Skip comments / blank lines
  if [[ "$file_path" =~ ^# ]] || [[ -z "$file_path" ]]; then
    continue
  fi

  input_file_name="$file_path"
  input_file_type="$file_type"
  input_beam_mode="${beam_mode:-UNKNOWN}"

  output_file_name="${output_dir}/xsec-ana-$(basename "${input_file_name}")"

  echo "----------------------------------------"
  echo "File ${counter}/${total_files}"
  echo "Input:  ${input_file_name}"
  echo "Type:   ${input_file_type}"
  echo "Beam:   ${input_beam_mode}"
  echo "Output: ${output_file_name}"
  echo "----------------------------------------"

  date
  time "${process_ntuples_bin}" \
    "${input_file_name}" \
    "${input_file_type}" \
    "${selections}" \
    "${output_file_name}" \
    "${input_beam_mode}"
  date

  # Write output file to analysis list
  echo "${output_file_name} ${input_file_type} ${input_beam_mode}" >> "${analysis_list_file}"

  counter=$((counter + 1))

done < "${ntuple_list_file}"

echo ""
echo "Done."
echo "Analysis file list saved to:"
echo "  ${analysis_list_file}"
