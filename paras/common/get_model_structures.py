import os
from shutil import copy
from sys import argv

from paras.common.fasta import read_fasta


def get_pdb_models(fasta, pdb_dir, structure_dir):
    if not os.path.exists(structure_dir):
        os.mkdir(structure_dir)
    domain_to_sequence = read_fasta(fasta)
    for file_name in os.listdir(pdb_dir):
        domain = file_name.split('_model')[0]
        if domain in domain_to_sequence:
            file_path = os.path.join(pdb_dir, file_name)
            file_destination = os.path.join(structure_dir, file_name)
            copy(file_path, file_destination)





if __name__ == "__main__":
    fasta = argv[1]
    pdb_dir = argv[2]
    out_dir = argv[3]

    get_pdb_models(fasta, pdb_dir, out_dir)


