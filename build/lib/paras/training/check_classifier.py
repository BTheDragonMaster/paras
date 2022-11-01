import os

import sklearn
from joblib import load

from paras.training.train_paras import train_paras
import paras.data

import paras.data.train_and_test_data.sequences
import paras.data.train_and_test_data.structure_alignments

TRAINING_SEQUENCES = os.path.join(os.path.dirname(paras.data.train_and_test_data.sequences.__file__),
                                  'train_sequences.fasta')
ALIGNMENT = os.path.join(os.path.dirname(paras.data.train_and_test_data.structure_alignments.__file__),
                         'training_alignment.fasta')

CLASSIFIER_DIR = os.path.dirname(paras.data.__file__)

CLASSIFIER = os.path.join(os.path.dirname(paras.data.__file__), 'paras.classifier')
CLASSIFIER_VERSION_DIR = os.path.join(os.path.dirname(paras.data.__file__), 'classifier_version.txt')


def get_classifier():
    with open(CLASSIFIER_VERSION_DIR, 'r') as classifier_version_file:
        classifier_version = classifier_version_file.read()

    scikit_learn_version = sklearn.__version__

    if not os.path.isfile(CLASSIFIER):
        print("Retraining PARAS classifier..")
        train_paras(TRAINING_SEQUENCES, ALIGNMENT, CLASSIFIER_DIR)
    elif classifier_version != scikit_learn_version:
        print("Retraining PARAS classifier..")
        train_paras(TRAINING_SEQUENCES, ALIGNMENT, CLASSIFIER_DIR)

    classifier = load(CLASSIFIER)
    with open(CLASSIFIER_VERSION_DIR, 'w') as classifier_version_file:
        classifier_version_file.write(scikit_learn_version)

    return classifier

