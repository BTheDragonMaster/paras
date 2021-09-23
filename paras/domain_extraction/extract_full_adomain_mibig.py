#!/usr/bin/env python

import subprocess
import os
from sys import argv
from pathlib import Path

from Bio import SearchIO

from paras.common.fasta import read_fasta

path = Path(os.getcwd())

HMM_DIR = f'{path}/data/AMP-binding_full.hmm'


def run_hmmscan(hmm_dir, fasta_dir, out_dir):
    program = 'hmmscan'
    arguments = [hmm_dir, fasta_dir]
    out_file = open(out_dir, 'w')

    command = [program] + arguments
    subprocess.call(command, stdout=out_file)
    out_file.close()


def remove_insertions(seq):
    new_seq = []
    for character in seq:
        if not character.islower():
            new_seq.append(character)

    new_seq = ''.join(new_seq)
    return new_seq


def make_header(ID, hit_id, start, end):
    return f'{ID}|{hit_id}|{start}-{end}'


def remove_gaps(sequence):
    new_sequence = []
    for character in sequence:
        if character != '-':
            new_sequence.append(character.upper())
    return ''.join(new_sequence)


def parse_hmm_results(hmm_results, fasta_out):
    fasta_file = open(fasta_out, 'w')
    for result in SearchIO.parse(hmm_results, 'hmmer3-text'):
        for hsp in result.hsps:
            if hsp.evalue < 0.00001:
                if hsp.hit_id == 'AMP-binding' or hsp.hit_id == 'AMP-binding_C':
                    seq = remove_gaps(hsp.query.seq)
                    header = make_header(result.id, hsp.hit_id, hsp.query_start, hsp.query_end)
                    fasta_file.write(">%s\n%s\n" % (header, seq))

    fasta_file.close()


def parse_fasta_id(fasta_id):
    id = fasta_id.split('|')[4]
    hit_id, hit_location = fasta_id.split('|')[-2:]
    print(id, hit_id, hit_location)
    hit_start, hit_end = hit_location.split('-')
    hit_start = int(hit_start)
    hit_end = int(hit_end)
    return id, hit_id, hit_start, hit_end


def find_amp_n_c_pairs(fasta_dir, original_fasta_dir, new_fasta_dir):
    fasta = read_fasta(fasta_dir)
    hits_by_seq_id = {}
    for id in fasta:
        seq_id, hit_id, hit_start, hit_end = parse_fasta_id(id)
        if not seq_id in hits_by_seq_id:
            hits_by_seq_id[seq_id] = []

        hits_by_seq_id[seq_id].append((hit_id, hit_start, hit_end))

    amp_nc_pairs = []

    print(hits_by_seq_id)

    for seq_id, hits in hits_by_seq_id.items():
        for hit_1 in hits:
            if hit_1[0] == 'AMP-binding':
                match_found = False
                for hit_2 in hits:
                    if hit_2[0] == 'AMP-binding_C':
                        if (hit_2[1] > hit_1[2]) and (hit_2[1] - hit_1[2] < 30):
                            match_found = True
                            amp_nc_pairs.append((seq_id, hit_1[1], hit_2[2]))

                if not match_found:
                    amp_nc_pairs.append((seq_id, hit_1[1], hit_1[2]))



    original_fasta = read_fasta(original_fasta_dir)
    amp_nc_pairs = sorted(amp_nc_pairs, key=lambda x: x[1])

    with open(new_fasta_dir, 'w') as new_fasta:
        for seq_id, sequence in original_fasta.items():
            counter = 0
            amp_id = seq_id.split('|')[4]
            for amp_nc_pair in amp_nc_pairs:
                pair_id, start, end = amp_nc_pair

                if amp_id == pair_id:
                    counter += 1

                    amp_sequence = sequence[start:end]
                    header = f'{amp_id}|{counter}|{start + 1}-{end}'
                    if amp_sequence:

                        new_fasta.write(f'>{header}\n{amp_sequence}\n')


if __name__ == "__main__":
    fasta = argv[1]
    file_label = '.'.join(fasta.split('.')[:-1])
    hmm_out = file_label + '.hmm_result'
    fasta_out = file_label + '_adomains_subdomains.faa'
    fasta_results = file_label + '_adomains.faa'

    run_hmmscan(HMM_DIR, fasta, hmm_out)
    parse_hmm_results(hmm_out, fasta_out)

    find_amp_n_c_pairs(fasta_out, fasta, fasta_results)


