from sys import argv

from paras.common.fasta import read_fasta, write_fasta

def rename_structural_alignment(fasta_file, out_file):
    domain_to_sequence = read_fasta(fasta_file)
    edited_domain_to_sequence = {}

    for domain, sequence in domain_to_sequence.items():
        domain_name = domain.split('_model')[0]
        edited_domain_to_sequence[domain_name] = sequence

    write_fasta(edited_domain_to_sequence, out_file)

def rename_trimmed_alignment(fasta_file, out_file):
    domain_to_sequence = read_fasta(fasta_file)
    edited_domain_to_sequence = {}

    for domain, sequence in domain_to_sequence.items():
        domain_name = domain.split(r'/')[0]
        edited_domain_to_sequence[domain_name] = sequence

    write_fasta(edited_domain_to_sequence, out_file)


if __name__ == "__main__":
    rename_structural_alignment(argv[1], argv[2])
