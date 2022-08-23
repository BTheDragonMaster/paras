from setuptools import setup, find_packages

VERSION = '0.0.1'
DESCRIPTION = 'PARAS: Predictive Algorithm for Resolving A-domain Specificity'
LONG_DESCRIPTION = 'Detect NRPS AMP-binding domains from an aa sequence and predict their substrate specificity'

setup(
    name="paras",
    version=VERSION,
    author="Barbara Terlouw",
    author_email="barbara.terlouw@wur.nl",
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    packages=find_packages(),
    package_data={"": ["AMP-binding_full.hmm*",
                       'Apositions*',
                       "*.fasta",
                       "*.txt",
                       "*.classifier"]},
    install_requires=[],
    scripts=['bin/paras', 'bin/paras-residues', 'bin/paras-test']
)
