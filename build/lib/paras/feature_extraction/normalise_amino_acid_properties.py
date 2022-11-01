from statistics import stdev, mean
from sys import argv

from paras.common.parsers import parse_unscaled_properties
from paras.common.writers import write_scaled_properties


def get_mean_and_stddev(aa_to_properties):
    std_devs = []
    means = []

    for i in range(len(aa_to_properties['A'])):
        column_values = [properties[i] for properties in aa_to_properties.values()]
        column_std_dev = stdev(column_values)
        column_mean = mean(column_values)
        std_devs.append(column_std_dev)
        means.append(column_mean)

    return std_devs, means


def normalise_aa_properties(properties_file, out_dir):
    aa_to_properties, categories = parse_unscaled_properties(properties_file)

    std_devs, means = get_mean_and_stddev(aa_to_properties)

    for aa, properties in aa_to_properties.items():

        for i, property in enumerate(properties):

            property_std_dev = std_devs[i]
            property_mean = means[i]

            aa_to_properties[aa][i] = (property - property_mean) / property_std_dev

    aa_to_properties['-'] = [0.0] * 15

    write_scaled_properties(aa_to_properties, categories, out_dir)


def get_mean_stddev(matrix):
    meanstddev = []

    for i in range(len(matrix[0])):
        col_vals = [row[i] for row in matrix]
        stddev = stdev(col_vals)
        mean_val = mean(col_vals)
        meanstddev.append((mean_val, stddev))

    return meanstddev


def normalise(matrix, meanstddev):
    for row in matrix:
        for i in range(len(row)):
            mean_val = meanstddev[i][0]
            stddev = meanstddev[i][1]
            row[i] = (row[i] - mean_val) / (stddev)


if __name__ == "__main__":
    normalise_aa_properties(argv[1], argv[2])
