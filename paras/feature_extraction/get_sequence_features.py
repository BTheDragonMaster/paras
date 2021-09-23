import os
import paras.data.aa_properties
from paras.common.parsers import parse_15_properties

PROPERTIES_FILE = os.path.join(os.path.dirname(paras.data.aa_properties.__file__), '15_aa_properties_normalised.txt')
PROPERTIES = parse_15_properties(PROPERTIES_FILE)


def get_sequence_features(sequence):
    features = []

    for aa in sequence:
        properties = PROPERTIES[aa]
        features += properties

    return features


def get_features(id_to_sequence):
    id_to_features = {}
    for seq_id, sequence in id_to_sequence.items():
        id_to_features[seq_id] = get_sequence_features(sequence)

    return id_to_features


def to_feature_vectors(id_to_features):
    ids = []
    feature_vectors = []
    for seq_id, features in id_to_features.items():
        ids.append(seq_id)
        feature_vectors.append(features)

    return ids, feature_vectors