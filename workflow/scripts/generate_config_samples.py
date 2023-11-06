import csv
import yaml
import os

# Get the current working directory
cwd = os.getcwd()

# Initialize the main dictionary to hold all species
main_dict = {'samples': {}}

# Read the CSV file and populate the dictionary
with open('../config/metadata_samples/samples.csv', 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        species = row['Species_long']
        sample_id = row['Sample_ID']
        # Check if species is already in the dictionary, if not, add it
        if species not in main_dict['samples']:
            main_dict['samples'][species] = {}
        # Now check if the sample_id is already in the species dict, if not, add it
        if sample_id not in main_dict['samples'][species]:
            main_dict['samples'][species][sample_id] = {
                'species': species,  # Include the species as an attribute
                'sample_id': sample_id,  # Include the species as an attribute
                'species_long': row['Species_long'],
                'sex': row['Sex'],
                'method': row['Method'],
                'orientation': row['Orientation'],
                'coverage': row['Coverage'],
                'path_prefix': row['Path_prefix'],
                'file_prefix': row['File_prefix'],
                'description': row['Description']
            }

# Output to YAML format
with open('../config/config_samples.yaml', 'w') as f:
    yaml.dump(main_dict, f, default_flow_style=False)
