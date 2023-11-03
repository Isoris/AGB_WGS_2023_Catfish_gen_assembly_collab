import csv
import yaml
import os

cwd = os.getcwd()

# Read the CSV file and populate a dictionary
samples = {}
with open('config/samples.csv', 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        sample_id = row['Sample_ID']
        species_short = row['Species']
        if species_short not in samples:
            samples[species_short] = {
                'sample_ids': [],
                'sexes': [],
                'methods': [],
                'orientations': [],
                'coverages': [],
                'path_prefixes': [],
                'file_prefixes': [],
                'species_long': [],
                'sex_long': [],
                'descriptions': []
            }
        samples[species_short]['sample_ids'].append(sample_id)
        samples[species_short]['sexes'].append(row['Sex'])
        samples[species_short]['methods'].append(row['Method'])
        samples[species_short]['orientations'].append(row['Orientation'])
        samples[species_short]['coverages'].append(row['Coverage'])
        samples[species_short]['path_prefixes'].append(row['Path_prefix'])
        samples[species_short]['file_prefixes'].append(row['File_prefix'])
        samples[species_short]['species_long'].append(row['Species_long'])
        samples[species_short]['sex_long'].append(row['Sex_long'])
        samples[species_short]['descriptions'].append(row['Description'])

# Output to YAML format
with open('config/config_samples.yaml', 'w') as f:
    yaml.dump({'samples': samples}, f, default_flow_style=False)


