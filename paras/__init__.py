#!/usr/bin/env python

import argparse

from paras.domain_extraction.extract_adomains import find_adomains

def define_arguments():
    parser = argparse.ArgumentParser(description="Run PARAS: Predictive Algorithm for Resolving A-domain Specificity")

    parser.add_argument("-f", "--fasta", type=str, required=True, help="Fasta file containing protein sequences")
    parser.add_argument("-o", "--out", type=str, required=True, help="Name of output directory")
    parser.add_argument("-g", "--group", action="store_true", help="Make group predictions in addition to single amino acid predictions")
    parser.add_argument("-s", "--stach", action="store_true", help="Save stachelhaus specificity conferring residues.")
    return parser

def run_paras(args):
    adomains_out = find_adomains(args.fasta, args.out)
    if args.stach:
        pass



