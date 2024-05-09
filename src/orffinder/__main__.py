from .orffinder import *

if __name__ == "__main__":
    if len(sys.argv) < 4:
        exit("Usage: python3 orffinder.py reference_genomes.fasta query_sequences.fasta [best|all]")
    genome, query = sys.argv[1], sys.argv[2]
    runtypes = ["best", "all"]
    runtype = sys.argv[3]
    if runtype not in runtypes:
        exit("Error: Unrecognised run type: " + runtype)
    if runtype == "best":
        OrfFinder(genome, query).bestOrf()
    elif runtype == "all":
        OrfFinder(genome, query).allGoodOrfs()
