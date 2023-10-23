#!/bin/bash

# Source the config.sh to get variables
source /path/to/config.sh

echo "Step 1: Adapter removal"
mkdir -p ${path_prefix_step1_trimming}AdaptRemoval

for y in "${!shortNames[@]}"; do
  if [[ -n "${reads1[y]}" && -n "${reads2[y]}" ]]; then
    ${path_prefix_programs}${progAdapterRemoval} \
      --file1 "${path_prefix_reads}${reads1[y]}" \
      --file2 "${path_prefix_reads}${reads2[y]}" \
      --gzip \
      --output1 "${path_prefix_reads_trimmed}${reads_trim1[y]}" \
      --output2 "${path_prefix_reads_trimmed}${reads_trim2[y]}" \
      --threads $CPUS
  fi
done

# Further steps would go here, each using their own respective prefixes for path locations
