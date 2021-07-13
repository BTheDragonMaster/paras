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


