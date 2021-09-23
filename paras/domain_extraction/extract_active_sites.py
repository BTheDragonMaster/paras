#!/usr/bin/env python

from paras.common.fasta import read_fasta, write_fasta

from typing import List, Tuple, Dict, Optional, Any, Union, IO
from pathlib import Path

import os
from sys import argv
from pprint import pprint
import subprocess

path = Path(os.getcwd())
parent = path.parent

ADOMAINS_FILENAME = f'{path}/data/A_domains_muscle.fasta'
APOSITION_FILENAME = f'{path}/data/Apositions_start_end.txt'
START_POSITION = 66

REF_SEQUENCE = "BAA00406.1.A1"


def run_muscle(in_file, out_file):

    command = ['muscle', '-quiet', '-profile',  '-in1', ADOMAINS_FILENAME, '-in2', f'{path}/{in_file}', '-out', f'{path}/{out_file}']

    subprocess.check_call(command)


def get_adomain_alignment(domain_name, domain_sequence):
    # Run muscle and collect sequence positions from file

    in_file = open('temp_in.common', 'w')

    in_file.write(f'>{domain_name}\n{domain_sequence}')
    in_file.close()

    run_muscle('temp_in.common', 'temp_out.common')
    alignments = read_fasta('temp_out.common')
    domain_alignment = alignments[domain_name]
    reference_alignment = alignments[REF_SEQUENCE]

    return reference_alignment, domain_alignment


#stolen from nrps_predictor.py
def read_positions(filename: str, start_position: int) -> List[int]:
    """ Loads positions from a tab-separated file. Positions are relative to the start_position.

        Arguments:
            filename: the path to the file containing the positions
            start_position: a relative start position to adjust all positions by

        Returns:
            a list of ints, one for each position found in the file
    """
    data = open(filename, "r")
    text = data.read().strip()
    results = []
    for i in text.split("\t"):
        results.append(int(i) - start_position)
    data.close()
    return results


def get_stach_aa_signature(reference_alignment: str, domain_alignment: str) -> str:
    """ Extract stachelhaus residues from A domains """
    positions = read_positions(APOSITION_FILENAME, START_POSITION)
    # Count residues in ref sequence and put positions in list
    poslist = build_position_list(positions, reference_alignment)
    # Extract positions from query sequence
    query_sig_seq = extract(domain_alignment, poslist)

    return query_sig_seq


#stolen from nrps_predictor.py
def build_position_list(positions: List[int], reference_seq: str) -> List[int]:
    """ Adjusts a list of positions to account for gaps in the reference sequence

        Arguments:
            positions: a list of ints that represent positions of interest in
                       the reference sequence
            reference_seq: the (aligned) reference sequence

        Returns:
            a new list of positions, each >= the original position
    """
    poslist = []
    position = 0
    for i, ref in enumerate(reference_seq):
        if ref != "-":
            if position in positions:
                poslist.append(i)
            position += 1
    return poslist


#stolen from nrps_predictor.py
def extract(sequence: str, positions: List[int]) -> str:
    """ Extracts a signature from an aligned sequence based on the provided
        positions. Accounts for gaps by looking behind or, if behind is already
        in the position list, ahead.

        Arguments:
            sequence: the aligned sequence to extract a signature from
            positions: the list of positions within the sequence to use

        Returns:
            the extracted signature as a string
    """

    start = positions[0]
    end = positions[1] + 1

    seq = sequence[start:end]

    gapless = []

    for aa in seq:
        if aa != '-':
            gapless.append(aa)

    return ''.join(gapless)


def write_tabular(id_to_signature, out_file):
    with open(out_file, 'w') as out:
        for id, signature in id_to_signature.items():
            out.write(f'{id}\t{signature}\n')


if __name__ == "__main__":
    fasta_file = argv[1]
    id_to_seq = read_fasta(fasta_file)
    id_to_signature = {}
    file_label = '.'.join(fasta_file.split('.')[:-1])
    out_table = file_label + '_stachelhaus.txt'

    for id, seq in id_to_seq.items():
        reference_alignment, domain_alignment = get_adomain_alignment(id, seq)
        signature = get_stach_aa_signature(reference_alignment, domain_alignment)
        print(id, signature)
        id_to_signature[id] = signature

    write_tabular(id_to_signature, out_table)
