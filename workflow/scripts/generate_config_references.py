import csv
import yaml

# Read the CSV file and populate a dictionary
references = {}
with open('../config/metadata_references/references.csv', 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        ref_id = row['Reference_ID']
        species_short = row['Species']
        if species_short not in references:
            references[species_short] = {
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
        references[species_short]['reference_ids'].append(ref_id)
        references[species_short]['species_long'].append(row['Species_long'])
        references[species_short]['sexes'].append(row['Sex'])
        references[species_short]['methods'].append(row['Method'])
        references[species_short]['orientations'].append(row['Orientation'])
        references[species_short]['descriptions'].append(row['Description'])
        references[species_short]['assembly_accessions'].append(row['Assembly_Accession'])
        references[species_short]['assembly_names'].append(row['Assembly_Name'])
        references[species_short]['organism_names'].append(row['Organism_Name'])
        for infraspecific_name in ['Breed', 'Strain', 'Cultivar', 'Ecotype', 'Isolate', 'Sex']:
            references[species_short]['organism_infraspecific_names'][infraspecific_name.lower()].append(row['Organism_Infraspecific_Names_' + infraspecific_name])
        references[species_short]['annotation_names'].append(row['Annotation_Name'])
        references[species_short]['assembly_stats_total_sequence_length'].append(row['Assembly_Stats_Total_Sequence_Length'])
        references[species_short]['assembly_level'].append(row['Assembly_Level'])
        references[species_short]['assembly_release_dates'].append(row['Assembly_Release_Date'])
        references[species_short]['wgs_project_accessions'].append(row['WGS_project_accession'])

# Output to YAML format
with open('../config/config_references.yaml', 'w') as f:
    yaml.dump({'references': references}, f, default_flow_style=False)

