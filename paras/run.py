
from typing import Tuple, List

from paras.domain_extraction.extract_residues import extract_stach_code
from paras.feature_extraction.get_sequence_features import get_sequence_features
from paras.training.check_classifier import get_classifier


def run_paras_bulk(sequences: str, threshold: float=0.5) -> List[List[Tuple[float, str]]]:
    classifier = get_classifier()
    amino_acid_classes = classifier.classes_
    all_probs_and_aas = []

    for sequence in sequences:

        stach, active_site = extract_stach_code(sequence)
        features = [get_sequence_features(active_site)]
        probabilities = classifier.predict_proba(features)[0]
        probs_and_aas = get_best_predictions(amino_acid_classes, probabilities, threshold)
        all_probs_and_aas.append(probs_and_aas)

    return all_probs_and_aas


def run_paras(sequence: str, threshold: float=0.5) -> List[Tuple[float, str]]:
    classifier = get_classifier()
    amino_acid_classes = classifier.classes_

    stach, active_site = extract_stach_code(sequence)
    features = [get_sequence_features(active_site)]
    probabilities = classifier.predict_proba(features)[0]
    probs_and_aas = get_best_predictions(amino_acid_classes, probabilities, threshold)

    return probs_and_aas


def get_best_predictions(amino_acid_classes, probabilities, threshold):
    probs_and_aa = []

    for i, probability in enumerate(probabilities):
        if probability >= threshold:
            substrate = amino_acid_classes[i].strip('"')
            probs_and_aa.append((probability, substrate))

    probs_and_aa.sort(reverse=True)

    return probs_and_aa


