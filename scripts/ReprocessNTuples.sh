#!/bin/bash

# Number of expected command-line arguments
num_expected=3

if [ "$#" -ne "$num_expected" ]; then
  echo "Usage: ./ReprocessNTuples.sh OUTPUT_DIRECTORY SELECTION_NAMES NTUPLE_LIST_FILE"
  exit 1
fi

output_dir=$1
selections=$2
ntuple_list_file=$3

if [ ! -f "$ntuple_list_file" ]; then
  echo "Ntuple list file \"${ntuple_list_file}\" not found"
  exit 1
fi

if [ ! -d "${output_dir}" ]; then
  echo "Output directory \"${output_dir}\" not found"
  exit 2
fi

total_files=$(grep -v '^\s*#' "$ntuple_list_file" | grep -v '^\s*$' | wc -l)
echo "Total number of files = ${total_files}"

counter=0
# Read file path and type from each non-comment line
while read -r file_path file_type; do
  # Skip comments and empty lines
  if [[ "$file_path" =~ ^# ]] || [[ -z "$file_path" ]]; then
    continue
  fi

  input_file_name="$file_path"
  input_file_type="$file_type"
  output_file_name="${output_dir}/xsec-ana-$(basename "${input_file_name}")"

  echo "Starting file: ${counter}/${total_files}"
  echo "Input file name: ${input_file_name}"
  echo "Input file type: ${input_file_type}"
  echo "Selections: ${selections}"
  echo "Output file name: ${output_file_name}"

  date
  time ProcessNTuples "${input_file_name}" "${input_file_type}" "${selections}" "${output_file_name}"
  date

  counter=$((counter + 1))
done < "${ntuple_list_file}"