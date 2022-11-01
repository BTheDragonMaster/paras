#!/usr/bin/env python

from typing import List

import os
from sys import argv

from paras.common.fasta import read_fasta, write_fasta
from paras.domain_extraction.read_positions import POSITIONS_STACH, POSITIONS_34, get_reference_positions
from paras.domain_extraction.align import align_adomain
import paras.data
import paras.data.train_and_test_data.structure_alignments

REF_SEQUENCE = "BAA00406.1.A1"
ADOMAIN_ALIGNMENT_FILE = os.path.join(os.path.dirname(paras.data.__file__), 'A_domains_structure_alignment.fasta')
ALIGNMENT_FILE = os.path.join(os.path.dirname(paras.data.train_and_test_data.structure_alignments.__file__), 'training_alignment.fasta')


def extract_stach_region(alignment, ref_sequence, out_file):
    id_to_seq = read_fasta(alignment)
    reference_alignment = id_to_seq[ref_sequence]
    position_start = POSITIONS_STACH[0]
    position_end = POSITIONS_STACH[-1]
    seq_range = [position_start, position_end]
    aligned_positions = get_reference_positions(seq_range, reference_alignment)

    id_to_truncated_seq = {}
    for id, seq in id_to_seq.items():
        id_to_truncated_seq[id] = seq[aligned_positions[0]:aligned_positions[1] + 1]

    write_fasta(remove_gaps(id_to_truncated_seq), out_file)


def stach_to_challis(id_to_stach):
    id_to_challis = {}
    for id, stach in id_to_stach.items():
        challis = stach[:8]
        id_to_challis[id] = challis

    return id_to_challis


def extract_stach_code(sequence):
    seq_id = "DOMAIN_TO_QUERY"
    aligned_id_to_seq = read_fasta(ALIGNMENT_FILE)
    if REF_SEQUENCE not in aligned_id_to_seq or seq_id not in aligned_id_to_seq:
        aligned_domain, aligned_reference = align_adomain(seq_id, sequence, ALIGNMENT_FILE)
    else:
        aligned_domain = aligned_id_to_seq[seq_id]
        aligned_reference = aligned_id_to_seq[REF_SEQUENCE]

    aligned_positions_stach = get_reference_positions(POSITIONS_STACH, aligned_reference)
    aligned_positions_34 = get_reference_positions(POSITIONS_34, aligned_reference)

    stach = []
    for position in aligned_positions_stach:
        stach.append(aligned_domain[position])
    stach = ''.join(stach)

    active_site = []
    for position in aligned_positions_34:
        active_site.append(aligned_domain[position])
    active_site = ''.join(active_site)
    return stach, active_site


def extract_stach_codes(fasta, alignment_file):

    id_to_seq = read_fasta(fasta)
    aligned_id_to_seq = read_fasta(alignment_file)

    id_to_stach = {}
    id_to_34 = {}

    for seq_id, seq in id_to_seq.items():
        if REF_SEQUENCE not in aligned_id_to_seq or seq_id not in aligned_id_to_seq:
            aligned_domain, aligned_reference = align_adomain(seq_id, seq, alignment_file)
        else:
            aligned_domain = aligned_id_to_seq[seq_id]
            aligned_reference = aligned_id_to_seq[REF_SEQUENCE]

        aligned_positions_stach = get_reference_positions(POSITIONS_STACH, aligned_reference)
        aligned_positions_34 = get_reference_positions(POSITIONS_34, aligned_reference)

        stach = []
        for position in aligned_positions_stach:
            stach.append(aligned_domain[position])
        stach = ''.join(stach)
        id_to_stach[seq_id] = stach

        active_site = []
        for position in aligned_positions_34:
            active_site.append(aligned_domain[position])
        active_site = ''.join(active_site)
        id_to_34[seq_id] = active_site

    return id_to_stach, id_to_34


def remove_gaps(id_to_seq):

    positions_to_remove = []

    alignment_length = len(id_to_seq[list(id_to_seq.keys())[0]])

    for position in range(alignment_length):
        only_gaps = True
        for seq_id, seq in id_to_seq.items():
            if seq[position] != '-':
                only_gaps = False
                break

        if only_gaps:
            positions_to_remove.append(position)

    cleaned_id_to_seq = {}

    for seq_id, seq in id_to_seq.items():
        cleaned_seq = []
        for position in range(alignment_length):
            if not position in positions_to_remove:
                cleaned_seq.append(seq[position])
        cleaned_seq = ''.join(cleaned_seq)
        cleaned_id_to_seq[seq_id] = cleaned_seq

    return cleaned_id_to_seq


def get_stach_aa_signature(reference_alignment: str, domain_alignment: str) -> str:
    """ Extract stachelhaus residues from A domains """

    # Count residues in ref sequence and put positions in list
    poslist = get_reference_positions(POSITIONS_STACH, reference_alignment)
    # Extract positions from query sequence
    query_sig_seq = extract(domain_alignment, poslist)
    # Add fixed lysine 517
    query_sig_seq += "K"

    return query_sig_seq


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
    seq = []
    for position in positions:
        aa = sequence[position]
        if aa == "-":
            if position - 1 not in positions:
                aa = sequence[position - 1]
            elif position + 1 not in positions:
                aa = sequence[position + 1]
        seq.append(aa)
    return "".join(seq)


def write_tabular(id_to_signature, out_file):
    with open(out_file, 'w') as out:
        for id, signature in id_to_signature.items():
            out.write(f'{id}\t{signature}\n')


if __name__ == "__main__":
    fasta = argv[1]
    alignment = argv[2]
    out = argv[3]
    id_to_stach, id_to_34 = extract_stach_codes(fasta, alignment)

    if not os.path.exists(out):
        os.mkdir(out)

    write_fasta(id_to_stach, os.path.join(out, 'stachelhaus.fasta'))
    write_fasta(id_to_34, os.path.join(out, 'active_site_34.fasta'))




