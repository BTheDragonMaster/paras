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
    if 'BAA00406.1.A1' in domain_to_specificity:
        del domain_to_specificity['BAA00406.1.A1']
    specificity_to_domain = reverse_dictionary(domain_to_specificity)
    train_domains = []
    test_domains = []

    for specificity, domains in specificity_to_domain.items():

        k = int(len(domains) * (1 - test_fraction))

        shuffle(domains)

        train_domains += domains[:k]
        test_domains += domains[k:]

    return train_domains, test_domains


def stratify_dataset_validation(fasta_file, fold, specificities_file, train_domains, out_dir):
    domain_to_specificity = parse_specificities(specificities_file)

    specificity_to_domain = {}
    for domain in train_domains:
        specificity = domain_to_specificity[domain]
        if specificity not in specificity_to_domain:
            specificity_to_domain[specificity] = []
        specificity_to_domain[specificity].append(domain)

    validation_training_sets = {}
    validation_test_sets = {}

    for i in range(fold):
        validation_training_sets[i] = []
        validation_test_sets[i] = []

    for specificity, domains in specificity_to_domain.items():
        shuffle(domains)

        for i, domain in enumerate(domains):
            validation_test_sets[i % fold].append(domain)
            for j in range(fold):
                if i % fold != j % fold:
                    validation_training_sets[j % fold].append(domain)

    domain_to_sequence = read_fasta(fasta_file)
    for i, training_domains in validation_training_sets.items():
        test_domains = validation_test_sets[i]
        validation_index = str(i + 1).zfill(2)

        training_domain_to_sequence = {}
        test_domain_to_sequence = {}

        for domain in training_domains:
            training_domain_to_sequence[domain] = domain_to_sequence[domain]

        for domain in test_domains:
            test_domain_to_sequence[domain] = domain_to_sequence[domain]

        out_directory = os.path.join(out_dir, f'validation_{validation_index}')

        out_training = os.path.join(out_directory, f'train_sequences.fasta')
        out_test = os.path.join(out_directory, f'validation_sequences.fasta')

        write_fasta(training_domain_to_sequence, out_training)
        write_fasta(test_domain_to_sequence, out_test)


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

    return train_domains

if __name__ == "__main__":
    sequences = argv[1]
    specificities = argv[2]
    test_fraction = float(argv[3])
    out_dir = argv[4]
    fold = int(argv[5])

    train_domains = split_train_test(sequences, specificities, test_fraction, out_dir)
    stratify_dataset_validation(sequences, fold, specificities, train_domains, out_dir)
