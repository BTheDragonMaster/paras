#!/usr/bin/env python

import os
from sys import argv
from Bio import SearchIO

from paras.domain_extraction.extract_adomains import run_hmmscan, HMM_FILE, TEMP_DIR
from paras.common.fasta import write_fasta
from paras.domain_extraction.read_positions import APOSITION_FILE_HMM, APOSITION_FILE_34_HMM, read_positions, get_reference_positions_hmm


def parse_hmm_results(hmm_results, out_dir):
    """
    Write sequences that match to AMP-binding domains to .fasta file

    Input:
    hmm_results: str, file location containing results of hmmscan
    fasta_out: str, file location for writing AMP-binding domains
    """

    domain_to_active_site = {}
    domain_to_stach = {}
    domain_to_score = {}
    domain_to_range = {}

    for result in SearchIO.parse(hmm_results, 'hmmer3-text'):
        for hsp in result.hsps:
            if hsp.evalue < 0.00001:
                if hsp.hit_id == 'AMP-binding':


                    if not hsp.query.id in domain_to_range:
                        domain_to_range[hsp.query.id] = [hsp.hit_start, hsp.hit_end]

                    if hsp.query.id not in domain_to_score or hsp.bitscore > domain_to_score[hsp.query.id]:
                        domain_to_score[hsp.query.id] = hsp.bitscore

                        stach_positions = get_reference_positions_hmm(read_positions(APOSITION_FILE_HMM, hsp.hit_start), hsp.hit.seq)
                        positions_34 = get_reference_positions_hmm(read_positions(APOSITION_FILE_34_HMM, hsp.hit_start), hsp.hit.seq)
                        stachelhaus = []
                        active_site = []

                        for position in stach_positions:
                            stachelhaus.append(hsp.query.seq[position])

                        for position in positions_34:
                            active_site.append(hsp.query.seq[position])

                        domain_to_active_site[hsp.query.id] = ''.join(active_site)
                        domain_to_stach[hsp.query.id] = ''.join(stachelhaus)

                    elif hsp.hit_start >= domain_to_range[hsp.query.id][1]:
                        stach_positions = get_reference_positions_hmm(read_positions(APOSITION_FILE_HMM, hsp.hit_start),
                                                                      hsp.hit.seq)
                        positions_34 = get_reference_positions_hmm(read_positions(APOSITION_FILE_34_HMM, hsp.hit_start),
                                                                   hsp.hit.seq)

                        for position in stach_positions:
                            domain_to_stach[hsp.query.id] += hsp.query.seq[position]

                        for position in positions_34:
                            domain_to_active_site[hsp.query.id] += hsp.query.seq[position]

                        domain_to_range[hsp.query.id][1] = hsp.hit_end

                    elif hsp.hit_end <= domain_to_range[hsp.query.id][0]:

                        stach_positions = get_reference_positions_hmm(read_positions(APOSITION_FILE_HMM, hsp.hit_start),
                                                                      hsp.hit.seq)
                        positions_34 = get_reference_positions_hmm(read_positions(APOSITION_FILE_34_HMM, hsp.hit_start),
                                                                   hsp.hit.seq)

                        for position in stach_positions:
                            domain_to_stach[hsp.query.id] += hsp.query.seq[position]

                        for position in positions_34:
                            domain_to_active_site[hsp.query.id] += hsp.query.seq[position]

                        stachelhaus = []
                        active_site = []

                        for position in stach_positions:
                            stachelhaus.append(hsp.query.seq[position])

                        for position in positions_34:
                            active_site.append(hsp.query.seq[position])

                        active_site = ''.join(active_site)
                        stachelhaus = ''.join(stachelhaus)

                        domain_to_active_site[hsp.query.id] = active_site + domain_to_active_site[hsp.query.id]
                        domain_to_stach[hsp.query.id] = stachelhaus + domain_to_stach[hsp.query.id]

                        domain_to_range[hsp.query.id][0] = hsp.hit_start

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    write_fasta(domain_to_active_site, os.path.join(out_dir, 'active_site_34.fasta'))
    write_fasta(domain_to_stach, os.path.join(out_dir, 'stachelhaus.fasta'))


if __name__ == "__main__":
    fasta = argv[1]
    out_dir = argv[2]
    hmm_out = os.path.join(TEMP_DIR, 'hmm_result.result')

    run_hmmscan(HMM_FILE, fasta, hmm_out)
    parse_hmm_results(hmm_out, out_dir)


