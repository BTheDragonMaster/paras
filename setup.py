import subprocess

from setuptools import setup, find_packages
from setuptools.command.install import install

VERSION = '0.0.3'
DESCRIPTION = 'PARAS: Predictive Algorithm for Resolving A-domain Specificity'
LONG_DESCRIPTION = 'Detect NRPS AMP-binding domains from an aa sequence and predict their substrate specificity'


class TrainParas(install):
    def run(self):
        #open("/Users/barbara/dummydummy", "w").write(str(glob.glob("bin/*")))
        subprocess.check_call(["python", "-c", "import os"])
        install.run(self)

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
    install_requires=["scikit-learn",
                      "biopython",
                      "joblib"],
    scripts=['bin/paras', 'bin/paras-residues', 'bin/paras-test']
)

