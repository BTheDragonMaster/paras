from sys import argv
from joblib import dump
import os

from sklearn.ensemble import RandomForestClassifier

from paras.domain_extraction.extract_residues import extract_stach_codes
from paras.common.fasta import write_fasta
from paras.feature_extraction.get_sequence_features import get_features, to_feature_vectors
from paras.common.parsers import parse_specificities
import paras.data

SPECIFICITIES_FILE = os.path.join(os.path.dirname(paras.data.__file__), 'specificities.txt')
SPECIFICITIES = parse_specificities(SPECIFICITIES_FILE)
REFERENCE_SEQUENCE = os.path.join(os.path.dirname(paras.data.__file__), 'reference_sequence.fasta')


def train_paras(sequence_file, alignment_file, out):
    id_to_stach, id_to_34 = extract_stach_codes(sequence_file, alignment_file)
    id_to_features = get_features(id_to_34)
    ids, feature_vectors = to_feature_vectors(id_to_features)
    response_vector = []
    for seq_id in ids:
        response_vector.append(SPECIFICITIES[seq_id])

    paras = RandomForestClassifier(random_state=25)
    paras.fit(feature_vectors, response_vector)

    if not os.path.exists(out):
        os.mkdir(out)

    dump(paras, os.path.join(out, 'paras.classifier'))
    write_fasta(id_to_stach, os.path.join(out, 'stach.fasta'))
    write_fasta(id_to_34, os.path.join(out, 'active_site_34.fasta'))


if __name__ == "__main__":
    sequences = argv[1]
    alignment = argv[2]
    out_dir = argv[3]

    train_paras(sequences, alignment, out_dir)







