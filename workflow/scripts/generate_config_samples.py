import csv
import yaml
import os

# Get the current working directory
cwd = os.getcwd()

# Initialize the main dictionary with a 'species' key
samples_dict = {'species': {}}

# Read the CSV file and populate the dictionary
with open('../config/metadata_samples/samples.csv', 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        sample_id = row['Sample_ID']
        species_short = row['Species']
        if species_short not in samples_dict['species']:
            samples_dict['species'][species_short] = {
                'sample_ids': [],
                'species_long': [],
                'sex': [],
                'method': [],
                'orientation': [],
                'coverage': [],
                'path_prefixes': [],
                'file_prefixes': [],
                'descriptions': []
            }
        species_data = samples_dict['species'][species_short]
        species_data['sample_ids'].append(sample_id)
        species_data['species_long'].append(row['Species_long'])
        species_data['sex'].append(row['Sex'])
        species_data['method'].append(row['Method'])
        species_data['orientation'].append(row['Orientation'])
        species_data['coverage'].append(row['Coverage'])
        species_data['path_prefixes'].append(row['Path_prefix'])
        species_data['file_prefixes'].append(row['File_prefix'])
        species_data['descriptions'].append(row['Description'])

# Output to YAML format
with open('../config/config_samples.yaml', 'w') as f:
    yaml.dump({'samples': samples_dict['species']}, f, default_flow_style=False)
