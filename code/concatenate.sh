#!/bin/bash

# Check if all four arguments are provided
if [ "$#" -lt 4 ]; then
  echo "Usage: $0 <input_string> <input_dir> <flag> <output_dir>"
  exit 1
fi

# Parse arguments
input_string="$1"
input_dir="$2"
base_folder="$3"
output_dir="$4"

# Ensure output directory exists
mkdir -p "$output_dir"

# Parse output file name and input sample names
IFS=":" read -r output_file input_files <<< "$input_string"
IFS="," read -ra input_file_array <<< "$input_files"

# Construct full output file name
output_file="$output_dir/$output_file.$base_folder.fastq.gz"

# Initialize array to hold valid fastq paths
fastq_files=()

# Locate matching input files
for value in "${input_file_array[@]}"; do
  for file in "$input_dir"/"$value".*.fastq.gz; do
    if [ -f "$file" ]; then
      fastq_files+=("$file")
    else
      echo "File not found: $file"
    fi
  done
done

# Fail if no valid input files
if [ ${#fastq_files[@]} -eq 0 ]; then
  echo "No valid FASTQ files found."
  exit 1
fi

# Concatenate with decompression and recompression
zcat "${fastq_files[@]}" | gzip > "$output_file"
echo "Concatenated ${#fastq_files[@]} files into $output_file"
