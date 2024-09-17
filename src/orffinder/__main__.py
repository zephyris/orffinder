from contextlib import nullcontext
from orffinder.gffoutput import GffOutput
from .orffinder import *
import argparse

def main():
    parser = argparse.ArgumentParser(description="Find open reading frames (ORFs) in the given sequence.")
    parser.add_argument('reference_genome', type=str, help="""The path to a FASTA file containing the reference genome used to cross-reference found ORFs against using tBLASTn.""")
    parser.add_argument('query_sequences', type=str, help="""The path to a FASTA file containing sequences ORFs should be found in. This could be e.g. a full genome or a set of mRNA sequences.""")
    parser.add_argument('--mode', type=str, choices=['best', 'all'], default='all', help="The mode of operation: 'best' finds the best ORF in a given sequence, while 'all' finds all good non-overlapping ORFs in a given sequence (see the documentation for details). Default is 'all'.")
    parser.add_argument('--output', type=str, default='-')

    args = parser.parse_args()
    orffinder = OrfFinder(args.reference_genome, args.query_sequences)
    if args.mode == "best":
        results = orffinder.bestOrf()
    elif args.mode == "all":
        results = orffinder.allGoodOrfs()

    gff = GffOutput(args.query_sequences)
    for orfs in results.values():
        for orf in orfs:
            gff.add_orf(orf)
    
    gff_str = gff.save()
    with open(args.output, 'w') if args.output != '-' else nullcontext(sys.stdout) as f:
        f.write(gff_str)

if __name__ == "__main__":
    main()
