# OrfFinder

A simple Python and tBLASTn-based open reading frame (ORF) finding tool.

## Installation

### Dependencies

OrfFinder requires `python`, and is easiest to install using `pip` and `git`. It also requires `ncbi-blast+` to be installed and available to the system command line.

On Linux, install using your normal package manager, for example:
``` shell
sudo apt update
sudo apt install ncbi-blast+ python3 python3-pip git
```

On Windows, download and install `ncbi-blast+` using the installer, available at [NCBI website](https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/).

You can install `python`, `pip` and `git` using `winget` in either Command Prompt or PowerShell:
``` shell
winget install Python.Python.3.0; winget install Git.Git
python3 -m pip install -update pip
```

### OrfFinder

Install using `pip` and `git`:
``` shell
pip install git+https://github.com/zephyris/orffinder
```

To reinstall and upgrade use `pip` and `git`:
``` shell
pip install --upgrade --force-reinstall git+https://github.com/zephyris/orffinder
```

To uninstall use `pip`
``` shell
pip uninstall orffinder
```

## Standalone usage

`orffinder` uses a few simple heuristics for finding ORFs.
It starts by finding stop codons, then extends them leftwards finding all possible start codons in frame with that stop codon, until it hits an in-frame stop codon or runs out of sequence.
These ORFs are then evaluated using simple properties (length, overlapping other ORFs, position) and tBLASTn against a reference genome collection.

Correct selection of start codon is then evaluated by tBLASTn of each start to start codon fragment of the ORF, extending short segments rightwards to a minimum length.
Start to start fragments of the ORFs with abnormally low tBLASTn hits for that protein are trimmed from the start of the ORF.

`orffinder` has two standard modes of operation:

1. Finding the best ORF in the forward direction in short sequence, like a transcript.
``` shell
python -m orffinder <reference_genomes.fasta> <query_transcript.fasta> best
```
This searches for the longest ORF, the ORF starting closest to the sequence start, and the ORF with the highest tBLASTn score against the reference genome collection. Then, start codons refined by tBLASTn.

2. Finding all good ORFs in either forward or reverse in a long sequence, like a genome or chromosome.
``` shell
python -m orffinder <reference_genomes.fasta> <query_genome.fasta> all
```
This searches for all ORFs, removing short ORFs which overlap a longer ORF, removing ORFs with low tBLASTn score against the reference genome collection. Then, start codons refined by tBLASTn.

## Python module usage

You can use `orffinder` in your Python scripts - however it is subject to change. The `OrfFinder` class carries out high-level operation. `Fasta` is used for data ingestion. `DnaSequence` and `Orf` are used for ORF finding and filtering. `Blast` handles tBLASTn searches.
