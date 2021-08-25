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

