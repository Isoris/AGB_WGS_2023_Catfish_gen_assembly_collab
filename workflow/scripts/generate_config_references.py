import csv
import yaml
import os
import re


# Read the existing YAML file to determine the starting point for ref_id_counter
def get_max_sample_id(yaml_file_path):
    with open(yaml_file_path, 'r') as f:
        config_samples = yaml.safe_load(f)
    
    max_id = 0
    for species, samples in config_samples['samples'].items():
        for sample_id in samples:
            # Assuming sample_id format is "ID<number>"
            current_id = int(re.search(r'ID(\d+)', sample_id).group(1))
            max_id = max(max_id, current_id)
    
    return max_id

# Initialize the ref_id_counter with the max ID from config_samples.yaml
ref_id_counter = get_max_sample_id('../config/config_samples.yaml')

# Continue with reading the references.csv and populating the references dictionary
references = {}

with open('../config/metadata_references/references.csv', 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        species = row['Species_long']
        if species not in references:
            references[species] = {
                # ... (rest of the initialization remains unchanged)
            }
        
        # Increment the counter to generate a unique ref_id for each entry
        ref_id_counter += 1
        ref_id = f"ID{ref_id_counter}"
        references[species]['reference_ids'].append(ref_id)

with open('../config/metadata_references/combined_references.csv', 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        species = row['Species_long']
        if species not in references:
            references[species] = {
                'reference_ids': [],
                'species_long': [],
                'sexes': [],
                'methods': [],
                'orientations': [],
                'descriptions': [],
                'assembly_accessions': [],
                'assembly_names': [],
                'organism_names': [],
                'organism_infraspecific_names': {
                    'breed': [],
                    'strain': [],
                    'cultivar': [],
                    'ecotype': [],
                    'isolate': [],
                    'sex': []
                },
                'annotation_names': [],
                'assembly_stats_total_sequence_length': [],
                'assembly_level': [],
                'assembly_release_dates': [],
                'wgs_project_accessions': []
            }
        
        # Increment the counter to generate a unique ref_id for each entry
        ref_id_counter += 1
        ref_id = f"ID{ref_id_counter}"
        references[species]['reference_ids'].append(ref_id)
        # ... (rest of your code remains unchanged)

# Output to YAML format
with open('../config/config_references.yaml', 'w') as f:
    yaml.dump({'references': references}, f, default_flow_style=False)
