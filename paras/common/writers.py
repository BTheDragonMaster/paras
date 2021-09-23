def write_tabular(dictionaries, labels, out_file):
    keys = sorted(list(dictionaries[0].keys()))
    header = 'Domain' + '\t'.join(labels) + '\n'

    with open(out_file, 'w') as out:
        out.write(header)
        for key in keys:
            out.write(key)
            for i, dictionary in enumerate(dictionaries):
                out.write(f'\t{dictionary[key]}')
            out.write('\n')


def write_scaled_properties(aa_to_properties, categories, out_file):
    aas = sorted(aa_to_properties.keys())
    with open(out_file, 'w') as out:
        out.write('AA')
        for category in categories:
            out.write(f'\t{category}')
        out.write('\n')
        for aa in aas:
            out.write(aa)
            properties = aa_to_properties[aa]

            for property in properties:
                out.write('\t%.10f' % property)

            out.write('\n')

