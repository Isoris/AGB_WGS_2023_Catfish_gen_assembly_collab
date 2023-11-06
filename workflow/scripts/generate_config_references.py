
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
        species = row.get('Species_long', 'Unknown Species')
        
        if species not in references:
            references[species] = {
                'reference_ids': [],
                'species_long': [],
                'sex': [],
                'method': [],
                'orientation': [],
                'common_names': [],
                'organism_infraspecific_names': {
                    'breed': [],
                    'strain': [],
                    'cultivar': [],
                    'ecotype': [],
                    'isolate': [],
                    'sex': []
                },                
                'qualifiers': [],
                'taxonomy_ids': [],
                'assembly_name': [],
                'assembly_accession': [],
                'sources': [],
                'annotations': [],
                'assembly_levels': [],
                'contig_n50s': [],
                'sizes': [],
                'submission_dates': [],
                'gene_counts': [],
                'bioprojects': [],
                'biosamples': [],
                'descriptions': [],

            }
        
        # Generate a unique ref_id for each entry and append
        ref_id_counter += 1
        ref_id = f"ID{ref_id_counter}"
        references[species]['reference_ids'].append(ref_id)

        # Append data using row.get for safety against missing columns
        references[species]['species_long'].append(row.get('Organism Scientific Name', ''))
        references[species]['common_names'].append(row.get('Organism Common Name', ''))
        references[species]['sex'].append(row.get('Sex', ''))
        references[species]['method'].append(row.get('Method', ''))
        references[species]['orientation'].append(row.get('Orientation', ''))
        
        # ... [do this for all infraspecific names and other attributes]
        
        references[species]['qualifiers'].append(row.get('Organism Qualifier', ''))
        references[species]['taxonomy_ids'].append(row.get('Taxonomy id', ''))
        references[species]['assembly_names'].append(row.get('Assembly Name', ''))
        references[species]['assembly_accessions'].append(row.get('Assembly Accession', ''))
        references[species]['sources'].append(row.get('Source', ''))
        references[species]['annotations'].append(row.get('Annotation', ''))
        references[species]['assembly_levels'].append(row.get('Level', ''))
        references[species]['contig_n50s'].append(row.get('Contig N50', ''))
        references[species]['sizes'].append(row.get('Size', ''))
        references[species]['submission_dates'].append(row.get('Submission Date', ''))
        references[species]['gene_counts'].append(row.get('Gene Count', ''))
        references[species]['bioprojects'].append(row.get('BioProject', ''))
        references[species]['biosamples'].append(row.get('BioSample', ''))
        references[species]['descriptions'].append(row.get('Description', ''))

# Output to YAML format
with open('../config/config_references.yaml', 'w') as f:
    yaml.dump({'references': references}, f, default_flow_style=False)