from typing import List
import os

import paras.data.positions

APOSITION_FILE = os.path.join(os.path.dirname(paras.data.positions.__file__), 'Apositions.txt')
APOSITION_FILE_34 = os.path.join(os.path.dirname(paras.data.positions.__file__), 'Apositions_34.txt')
APOSITION_FILE_HMM = os.path.join(os.path.dirname(paras.data.positions.__file__), 'Apositions_hmm.txt')
APOSITION_FILE_34_HMM = os.path.join(os.path.dirname(paras.data.positions.__file__), 'Apositions_34_hmm.txt')
START_POSITION = 66


def read_positions(position_file: str, start_position: int) -> List[int]:
    """
    Return positions from a tab-separated file. Positions are relative to start_position.

    Input:
    position_file: str, the path to tab-separated file containing the positions
    start_position: int, a relative start position to adjust all positions by

    Output:
    positions: list of [int, ->], one for each position found in the file
    """

    with open(position_file, 'r') as position_data:
        text = position_data.read().strip()
        positions = []
        for i in text.split("\t"):
            positions.append(int(i) - start_position)
    return positions


def get_reference_positions(positions: List[int], aligned_reference: str) -> List[int]:
    """
    Adjusts a list of positions to account for gaps in the reference sequence

    Input:
    positions: list of [int, ->], with integers representing positions of interest in
        the reference sequence
    aligned_reference: the (aligned) reference sequence

    Output:
    pos_list: list of [int, ->], a new list of positions, each >= the original position
    """
    pos_list = []
    position = 0
    for i, aa in enumerate(aligned_reference):
        if aa != "-":
            if position in positions:
                pos_list.append(i)
            position += 1
    return pos_list


def get_reference_positions_hmm(positions: List[int], aligned_reference: str) -> List[int]:
    """
    Adjusts a list of positions to account for gaps in the reference sequence

    Input:
    positions: list of [int, ->], with integers representing positions of interest in
        the reference sequence
    aligned_reference: the (aligned) reference sequence

    Output:
    pos_list: list of [int, ->], a new list of positions, each >= the original position
    """
    pos_list = []
    position = 0
    for i, aa in enumerate(aligned_reference):
        if aa != ".":
            if position in positions:
                pos_list.append(i)
            position += 1

    return pos_list


POSITIONS_STACH = read_positions(APOSITION_FILE, START_POSITION)
POSITIONS_34 = read_positions(APOSITION_FILE_34, START_POSITION)

