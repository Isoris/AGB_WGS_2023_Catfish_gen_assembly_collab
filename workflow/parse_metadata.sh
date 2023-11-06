#!/bin/bash

# Set the working directory to the script's directory
cd scripts/

# Activate the quality control conda environment if necessary
# source activate quality_control

# Define the output YAML files
config_samples="config_samples.yaml"
config_references="config_references.yaml"

# Define the log files
log_samples="generate_config_samples.log"
log_references="generate_config_references.log"

# Run Python scripts to generate config YAMLs and log the outputs
echo "Generating $config_samples..."
python generate_config_samples.py > "../config/$config_samples" 2> "$log_samples"

cd ../config/metadata_references

cat *.ref.tsv | sed 's/\t/,/g' > combined_references.csv

cd ../../scripts

echo "Generating $config_references..."
python generate_config_references.py > "../config/$config_references" 2> "$log_references"

# Concatenate the generated YAML files with the base config.yaml
echo "Merging configuration files into config.yaml..."
cat ../config/config_base.yaml ../config/${config_samples} ../config/${config_references} > ../config/config.yaml

# Output the resulting config.yaml to the console (optional)
echo "Merged config.yaml:"
head -n 10 ../config/config.yaml

# Return to the original directory
cd ..

