
##### Definition of functions

def get_unique_values(key):
    values = set()
    for species in config["samples"].values():
        for sample in species.values():
            value = sample.get(key)
            if value is not None:
                values.add(value)
    return list(values)

def get_unique_values_references(key):
    values = set()
    for sample in config["references"]["references"].values():
        values.update(sample.get(key, []))
    return list(values)

def get_sample_ids():
    return get_unique_values('sample_id')

def get_reference_accession():
    return get_unique_values('assembly_accession')

def get_species_file_prefix():
    return get_unique_values('file_prefixes')

def get_species_long():
    return get_unique_values('species_long')

def get_species():
    return get_unique_values('species')

def get_sex():
    return get_unique_values('sex')

def get_method():
    return get_unique_values('method')

def get_orientation():
    return get_unique_values('orientation')
