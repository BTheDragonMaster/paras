import os

from paras.domain_extraction.read_positions import POSITIONS_34, POSITIONS_STACH
from paras.common.writers import write_tabular
from paras.common.fasta import read_fasta

from sys import argv


def compare_sequences(sequence_1, sequence_2, type='stach'):
    unmatching_positions = []
    mismatches = []

    if type == 'stach':
        positions = POSITIONS_STACH
    elif type == '34':
        positions = POSITIONS_34

    for i, aa_1 in enumerate(sequence_1):
        try:
            aa_2 = sequence_2[i]

            if aa_1 != aa_2:
                position = positions[i]
                unmatching_positions.append(position)
                mismatches.append((aa_1, aa_2))

        except IndexError:
            pass

    return unmatching_positions, mismatches


def compare_alignment_methods(sequence_based_extraction, structure_based_extraction, extraction_type, out_dir):

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    domain_to_sequence = read_fasta(sequence_based_extraction)
    domain_to_structure = read_fasta(structure_based_extraction)

    domain_to_mismatches = {}

    domain_to_mismatching_positions = {}

    # Keeps track of how many gaps there are in total in all sequences
    gaps_sequence = 0
    gaps_structure = 0

    # Keeps track of how many sequences contain gaps
    gapped_sequence = 0
    gapped_structure = 0

    domain_to_gap_count_sequence = {}
    domain_to_gap_count_structure = {}

    for domain, sequence in domain_to_sequence.items():
        structure = domain_to_structure[domain]

        gap_count_sequence = sequence.count('-')
        gap_count_structure = structure.count('-')

        if gap_count_sequence:
            gapped_sequence += 1
            gaps_sequence += gap_count_sequence

        if gap_count_structure:
            gapped_structure += 1
            gaps_structure += gap_count_structure

        domain_to_gap_count_sequence[domain] = gap_count_sequence
        domain_to_gap_count_structure[domain] = gap_count_structure

        unmatching_positions, mismatches = compare_sequences(sequence, structure, type=extraction_type)

        domain_to_mismatching_positions[domain] = ', '.join(map(str, unmatching_positions))
        domain_to_mismatches[domain] = ', '.join([x[0] + x[1] for x in mismatches])

    write_tabular([domain_to_sequence,
                   domain_to_structure,
                   domain_to_gap_count_sequence,
                   domain_to_gap_count_structure,
                   domain_to_mismatching_positions,
                   domain_to_mismatches],
                  ['sequence_based_extraction',
                   'hmm_based_extraction',
                   'gap_count_sequence_based',
                   'gap_count_structure_based',
                   'mismatching_positions',
                   'mismatches'],
                  os.path.join(out_dir, f'{extraction_type}_comparison.txt'))


if __name__ == "__main__":
    sequence_based_extraction, structure_based_extraction = argv[1], argv[2]
    extraction_type = argv[3]
    out_dir = argv[4]

    compare_alignment_methods(sequence_based_extraction, structure_based_extraction, extraction_type, out_dir)


