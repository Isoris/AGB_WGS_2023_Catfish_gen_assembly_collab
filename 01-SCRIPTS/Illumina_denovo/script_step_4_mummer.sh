#!/bin/bash

# Source the config file
source /path/to/config.sh

# Set the current step for output directory
current_step="mummer"
output_directory=${output[${current_step}]}

echo "Starting to run MuMmer for comparing the ${Species}_${Sex}_${Version} or ${shortNames}"

# Ensure the output directory exists
mkdir -p ${output_directory}

# Loop through the reference genomes
for ref_genome in "${ref_genome_clarias_macrocephalus}" "${ref_genome_clarias_gariepinus}" "${ref_genome_f1}"; do
  # Check if reference genome exists
  if [[ -f ${ref_genome} ]]; then
    ref_genome_name=$(basename -- "${ref_genome}")
    ref_genome_name_no_ext="${ref_genome_name%.*}"
    
    # Run MUMmer to compare assemblies
    nucmer -p ${output_directory}/${ref_genome_name_no_ext}_output assembly_unicycler.fasta ${ref_genome}
    
    # Run dnadiff for additional statistics (optional)
    dnadiff -p ${output_directory}/${ref_genome_name_no_ext}_dnadiff ${ref_genome} assembly_unicycler.fasta
    
    # Generate a plot of the comparison
    mummerplot -p ${output_directory}/${ref_genome_name_no_ext}_output ${output_directory}/${ref_genome_name_no_ext}_output.delta
  else
    echo "Reference genome ${ref_genome} not found. Skipping..."
  fi
done

echo "MuMmer comparison completed"

