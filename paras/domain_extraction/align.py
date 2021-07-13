import os
import subprocess

import paras.data
import paras.data.temp
from paras.common.fasta import read_fasta

#ADOMAIN_ALIGNMENT_FILE = os.path.join(os.path.dirname(paras.data.__file__), 'A_domains_muscle.fasta')
TEMP_DIR = os.path.dirname(paras.data.temp.__file__)
#REF_SEQUENCE = "P0C062_A1"

REF_SEQUENCE = "BAA00406.1.A1"
ADOMAIN_STRUCTURAL_ALIGNMENT_FILE = os.path.join(os.path.dirname(paras.data.__file__), 'A_domains_structure_alignment.fasta')
ADOMAIN_SEQUENCE_ALIGNMENT_FILE = os.path.join(os.path.dirname(paras.data.__file__), 'A_domains_structure_alignment.fasta')


def run_muscle(in_file, alignment_file, out_file):

    command = ['muscle', '-quiet', '-profile',  '-in1', alignment_file, '-in2', in_file, '-out', out_file]

    subprocess.check_call(command)

def align_adomain(domain_name, domain_sequence, alignment_file):
    # Run muscle and collect sequence positions from file

    temp_in = os.path.join(TEMP_DIR, 'temp_in_alignment.fasta')
    temp_out = os.path.join(TEMP_DIR, 'temp_out_alignment.fasta')

    with open(temp_in, 'w') as temp:
        temp.write(f'>{domain_name}\n{domain_sequence}')

    #align sequence to all a domains in database with muscle
    run_muscle(temp_in, alignment_file, temp_out)
    id_to_alignment = read_fasta(temp_out)

    #aligned sequence of domain
    aligned_domain = id_to_alignment[domain_name]

    #aligned sequence of 1AMU reference sequence
    aligned_reference = id_to_alignment[REF_SEQUENCE]

    return aligned_domain, aligned_reference