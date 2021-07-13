#!/usr/bin/env python

import os
from sys import argv
from random import shuffle

from paras.common.parsers import parse_specificities
from paras.common.fasta import read_fasta, write_fasta

def reverse_dictionary(dictionary):
    reverse_dict = {}
    for key, value in dictionary.items():
        if value not in reverse_dict:
            reverse_dict[value] = []
        reverse_dict[value].append(key)

    return reverse_dict


def stratify_dataset(test_fraction, domain_to_specificity):
    specificity_to_domain = reverse_dictionary(domain_to_specificity)
    train_domains = []
    test_domains = []

    for specificity, domains in specificity_to_domain.items():

        k = int(len(domains) * (1 - test_fraction))

        shuffle(domains)

        train_domains += domains[:k]
        test_domains += domains[k:]

    return train_domains, test_domains


def split_train_test(fasta_file, specificities_file, test_fraction, out_dir):
    domain_to_specificity = parse_specificities(specificities_file)
    domain_to_sequence = read_fasta(fasta_file)

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    train_sequences = os.path.join(out_dir, 'train_sequences.fasta')
    test_sequences = os.path.join(out_dir, 'test_sequences.fasta')

    train_domains, test_domains = stratify_dataset(test_fraction, domain_to_specificity)

    train_domain_to_sequence = {}
    for domain in train_domains:
        train_domain_to_sequence[domain] = domain_to_sequence[domain]

    test_domain_to_sequence = {}

    for domain in test_domains:
        test_domain_to_sequence[domain] = domain_to_sequence[domain]

    write_fasta(train_domain_to_sequence, train_sequences)
    write_fasta(test_domain_to_sequence, test_sequences)

if __name__ == "__main__":
    sequences = argv[1]
    specificities = argv[2]
    test_fraction = float(argv[3])
    out_dir = argv[4]

    split_train_test(sequences, specificities, test_fraction, out_dir)
