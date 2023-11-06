
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
        references[species]['species_long'].append(row['species_long'])
        references[species]['sexes'].append(row['Sex'])
        references[species]['methods'].append(row['Method'])
        references[species]['orientations'].append(row['Orientation'])
        references[species]['descriptions'].append(row['Description'])
        references[species]['assembly_accessions'].append(row['Assembly_Accession'])
        references[species]['assembly_names'].append(row['Assembly_Name'])
        references[species]['organism_names'].append(row['Organism_Name'])
        for infraspecific_name in ['Breed', 'Strain', 'Cultivar', 'Ecotype', 'Isolate', 'Sex']:
            if row.get(f'Organism_Infraspecific_Names_{infraspecific_name}'):
                references[species]['organism_infraspecific_names'][infraspecific_name.lower()].append(
                    row[f'Organism_Infraspecific_Names_{infraspecific_name}']
                )
        references[species]['annotation_names'].append(row['Annotation_Name'])
        references[species]['assembly_stats_total_sequence_length'].append(row['Assembly_Stats_Total_Sequence_Length'])
        references[species]['assembly_level'].append(row['Assembly_Level'])
        references[species]['assembly_release_dates'].append(row['Assembly_Release_Date'])
        references[species]['wgs_project_accessions'].append(row['WGS_project_accession'])

# Output to YAML format
with open('../config/config_references.yaml', 'w') as f:
    yaml.dump({'references': references}, f, default_flow_style=False)
