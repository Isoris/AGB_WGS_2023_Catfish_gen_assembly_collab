#!/bin/bash

# Source the config.sh to get variables
source /path/to/config.sh

# Step 0: Quality control with FastQC
output_dir="${output[fastqc]}"
mkdir -p $output_dir

# Function to check if FastQC has already been run
check_fastqc_run () {
  local read_file="$1"
  local output_check="${output_dir}/${read_file}_fastqc.zip"

  if [[ -e "$output_check" ]]; then
    echo "FastQC already run on ${read_file}. Skipping..."
    return 1
  else
    return 0
  fi
}

for sample in "${!shortNames[@]}"; do
  for r in "reads1" "reads2"; do
    if [[ -n "${!r[sample]}" ]]; then
      read_file="${!r[sample]}"

      # Check if FastQC already run on raw reads
      if check_fastqc_run "$read_file"; then
        ${progfastqc} -o $output_dir -t $CPUS "${path_prefix_reads}${read_file}"
      fi

      # If trimmed reads exist, run FastQC on them too
      read_file_trimmed="${reads_trim1[sample]}"
      if [[ -e "${path_prefix_reads_trimmed}${read_file_trimmed}" ]]; then
        if check_fastqc_run "$read_file_trimmed"; then
          ${progfastqc} -o $output_dir -t $CPUS "${path_prefix_reads_trimmed}${read_file_trimmed}"
        fi
      fi
    fi
  done
done

echo "Step 0: FastQC completed."

