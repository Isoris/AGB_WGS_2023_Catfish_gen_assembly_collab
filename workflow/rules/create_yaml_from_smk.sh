#!/bin/bash

# Loop through all .smk files in the current directory
for smk_file in *.smk; do
  # Extract the base name of the .smk file to use for the .yaml file
  base_name=$(basename -- "$smk_file" .smk)
  yaml_file="${base_name}.yaml"
  
  # Create a new .yaml file if it doesn't already exist
  if [ ! -f "$yaml_file" ]; then
    echo "Creating $yaml_file..."
    
    # Add default content to the .yaml file
    echo "dependencies:" > "$yaml_file"
    echo "  - package1" >> "$yaml_file"
    echo "  - package2" >> "$yaml_file"
  fi
done
