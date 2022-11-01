def parse_smiles(smiles_file, lowercase=False):
    substrate_to_smiles = {}
    with open(smiles_file, 'r') as names_and_smiles:
        for line in names_and_smiles:
            line = line.strip()
            name, smiles = line.split()
            if lowercase:
                substrate_to_smiles[name.lower()] = smiles
            else:
                substrate_to_smiles[name] = smiles

    return substrate_to_smiles


def parse_specificities(specificities_file):
    domain_to_specificity = {}
    with open(specificities_file, 'r') as specificities:
        for line in specificities:
            line = line.strip()
            if line:
                domain, specificity = line.split('\t')
                specificity = specificity.lower()
                domain_to_specificity[domain] = specificity

    return domain_to_specificity


def parse_15_properties(properties_file):
    aa_to_vector = {}

    with open(properties_file, 'r') as properties:

        properties.readline()

        for line in properties:
            line = line.strip()
            aa = line.split()[0]
            properties = line.split()[1:]
            properties_float = []
            for property in properties:
                properties_float.append(float(property))
            assert len(properties_float) == 15
            aa_to_vector[aa] = properties_float

    return aa_to_vector


def parse_unscaled_properties(properties_file):
    aa_to_vector = {}

    with open(properties_file, 'r') as properties:
        header = properties.readline()
        categories = header.split('\t')[2:]
        categories_cleaned = []
        for category in categories:
            categories_cleaned.append(category.strip())

        for line in properties:
            line = line.strip()
            aa = line.split('\t')[1]
            properties = line.split('\t')[2:]
            properties_float = []
            for property in properties:
                properties_float.append(float(property))
            assert len(properties_float) == 15
            aa_to_vector[aa] = properties_float

    return aa_to_vector, categories_cleaned
