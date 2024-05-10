from .orffinder import *
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find open reading frames (ORFs) in the given sequence.")
    parser.add_argument('reference_genome', type=str, help="""The path to a FASTA file containing the reference genome used to cross-reference found ORFs against using tBLASTn.""")
    parser.add_argument('query_sequences', type=str, help="""The path to a FASTA file containing sequences ORFs should be found in. This could be e.g. a full genome or a set of mRNA sequences.""")
    parser.add_argument('--mode', type=str, choices=['best', 'all'], default='all', help="The mode of operation: 'best' finds the best ORF in a given sequence, while 'all' finds all possible non-overlapping ORFs in a given sequence. Default is 'all'.")

    args = parser.parse_args()
    print()
    orffinder = OrfFinder(args.reference_genome, args.query_sequences)
    if args.mode == "best":
        results = orffinder.bestOrf()
    elif args.mode == "all":
        results = orffinder.allGoodOrfs()
