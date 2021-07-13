def write_tabular(dictionaries, labels, out_file):
    keys = sorted(list(dictionaries[0].keys()))
    header = '\t'.join(labels) + '\n'

    with open(out_file, 'w') as out:
        out.write(header)
        for key in keys:
            for i, dictionary in enumerate(dictionaries):
                if i == 0:
                    out.write(dictionary[key])
                else:
                    out.write(f'\t{dictionary[key]}')
            out.write('\n')

