from paras.domain_extraction.read_positions import POSITIONS_34, POSITIONS_STACH
from paras.common.writers import write_tabular



def compare_sequences(sequence_1, sequence_2, type='stach'):
    unmatching_positions = []
    mismatches = []

    if type == 'stach':
        positions = POSITIONS_STACH
    elif type == '34':
        positions = POSITIONS_34

    for i, aa_1 in enumerate(sequence_1):
        aa_2 = sequence_2[i]

        if aa_1 != aa_2:
            position = positions[i]
            unmatching_positions.append(position)
            mismatches.append((aa_1, aa_2))

    return unmatching_positions, mismatches



