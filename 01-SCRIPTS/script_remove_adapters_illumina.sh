#!/bin/bash

readsPath="./" # Add your default path here
progAdapterRemoval="/path/to/AdapterRemoval" # Path to AdapterRemoval executable

# Initialize empty arrays
CPU
shortNames=()
reads1=()
reads2=()

# Read from the TSV file and populate arrays
while IFS=$'\t' read -r species_short method prefix orientation path files coverage_x sex species; do
  if [[ "$method" == "Illumina_PE" ]]; then
    key="${species_short}_${prefix}_${orientation}"
    fullPath="${path}${files}"

    shortNames+=("$key")

    if [[ "$orientation" == "R1" ]]; then
      reads1+=("$fullPath")
    elif [[ "$orientation" == "R2" ]]; then
      reads2+=("$fullPath")
    fi
  fi
done < <(tail -n +2 "read_data.tsv")  # Assuming your TSV is named data.tsv and skipping the header line

# Step 3: Adapter removal
echo "Step 3: Adapter removal"
mkdir -p AdaptRemoval

for y in "${!shortNames[@]}"; do
  if [[ -n "${reads1[y]}" && -n "${reads2[y]}" ]]; then  # Make sure both R1 and R2 are set
    ${progAdapterRemoval}/AdapterRemoval \
      --file1 "${reads1[y]}" \
      --file2 "${reads2[y]}" \
      --gzip \
      --output1 "AdaptRemoval/${shortNames[y]}_R1_TRIM.gz" \
      --output2 "AdaptRemoval/${shortNames[y]}_R2_TRIM.gz" \
      --threads 8
  fi
done

