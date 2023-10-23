import csv
import yaml

# Read the CSV file and populate a dictionary
samples = {}
with open('../../config/samples.csv', 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        sample_id = row['species_short']
        if sample_id not in samples:
            samples[sample_id] = {
                'methods': [],
                'prefixes': [],
                'orientations': [],
                'paths': [],
                'files': [],
                'descriptions': [],
                'coverages': [],
                'sex': row['sex'],
                'species': row['species'],
                'outputs': {}
            }
        samples[sample_id]['methods'].append(row['method'])
        samples[sample_id]['prefixes'].append(row['prefix'])
        samples[sample_id]['orientations'].append(row['orientation'])
        samples[sample_id]['paths'].append(row['path'])
        samples[sample_id]['files'].append(row['files'])
        samples[sample_id]['descriptions'].append(row['description'])
        samples[sample_id]['coverages'].append(row['coverage_X'])

# Output to YAML format
with open('../../config/config_samples.yaml', 'w') as f:
    yaml.dump({'samples': samples}, f, default_flow_style=False)


