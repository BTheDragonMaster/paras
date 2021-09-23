from sys import argv
from shutil import copy
import os

from paras.common.fasta import read_fasta


def copy_structures(structures_dir, sequence_dir, out_dir):
    train_to_sequence = read_fasta(os.path.join(sequence_dir, 'train_sequences.fasta'))

    train_structures_dir = os.path.join(out_dir, 'train')
    validation_structures_dir = os.path.join(out_dir, 'validation')


    if not os.path.exists(train_structures_dir):
        os.mkdir(train_structures_dir)

    if not os.path.exists(validation_structures_dir):
        os.mkdir(validation_structures_dir)

    reference_structure = os.path.join(structures_dir, 'BAA00406.1.A1_model.pdb')
    copy(reference_structure, os.path.join(train_structures_dir, 'BAA00406.1.A1_model.pdb'))

    for domain in train_to_sequence:
        domain_file = os.path.join(structures_dir, f'{domain}_model.pdb')
        if os.path.exists(domain_file):
            copy(domain_file, os.path.join(train_structures_dir, f'{domain}_model.pdb'))
        else:
            print(f"No model for {domain}")

    for folder in os.listdir(sequence_dir):
        if folder.startswith('validation'):
            validation_dir = os.path.join(sequence_dir, folder)

            validation_structure_dir = os.path.join(validation_structures_dir, folder)
            if not os.path.exists(validation_structure_dir):
                os.mkdir(validation_structure_dir)

            train_validation_to_sequence = read_fasta(os.path.join(validation_dir, 'train_sequences.fasta'))
            for domain in train_validation_to_sequence:
                structure = os.path.join(structures_dir, f'{domain}_model.pdb')
                if os.path.exists(structure):
                    copy(structure, os.path.join(validation_structure_dir, f'{domain}_model.pdb'))
                else:
                    print("Oops")
                

            copy(reference_structure, os.path.join(validation_structure_dir, 'BAA00406.1.A1_model.pdb'))


if __name__ == "__main__":
    structures = argv[1]
    sequences = argv[2]
    out = argv[3]

    copy_structures(structures, sequences, out)

