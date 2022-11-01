import os
from typing import Dict
import paras.data

SMILES_FILE = os.path.join(os.path.dirname(paras.data.__file__), 'smiles.txt')


def load_smiles() -> Dict[str, str]:
    name_to_smiles = {}
    with open(SMILES_FILE, 'r') as smiles_file:
        smiles_file.readline()
        for line in smiles_file:
            name, smiles = line.strip().split('\t')
            name_to_smiles[name] = smiles
    return name_to_smiles


UNKNOWN_SUBSTRATE = "**Unknown**"

_SMILES = load_smiles()
_SMILES[UNKNOWN_SUBSTRATE] = "NC(*)C(=O)O"


def get_smiles(substrate_name: str) -> str:
    """
    Return SMILES string from substrate name

    substrate_name: str. Unknown substrate: use UNKNOWN_SUBSTRATE
    """
    if not isinstance(substrate_name, str):
        raise TypeError(f"Expected str, got {type(substrate_name)}.")

    if substrate_name in _SMILES:
        return _SMILES[substrate_name]

    raise ValueError(f"Cannot fetch SMILES string for {substrate_name}.")
