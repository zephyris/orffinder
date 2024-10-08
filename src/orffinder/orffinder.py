import os, sys, subprocess, statistics, multiprocessing
from functools import cached_property
from typing import Union

class BlastHit:
    """
    Convenience function for handling BLAST -outfmt 6 data lines
    """
    def __init__(self, outfmt6_line):
        self.data = outfmt6_line.split("\t")
        self.qseqid = self.data[0]
        self.sseqid = self.data[1]
        self.pident = float(self.data[2])
        self.length = int(self.data[3])
        self.mismatch = int(self.data[4])
        self.gapopen = int(self.data[5])
        self.qstart = int(self.data[6])
        self.qend = int(self.data[7])
        self.sstart = int(self.data[8])
        self.send = int(self.data[9])
        self.evalue = float(self.data[10])
        self.bitscore = float(self.data[11])

    def __repr__(self):
        return "\t".join(self.data)

class Blast:
    """
    Object to handle BLAST indexing and searching.
    ncbi-blast+ must be installed and the various programs available in path.
    Query using lists of objects with .fasta() and .sequenceName() methods.
    """
    def __init__(self, subject_path, db_type, program_type, max_hits=None, do_indexing=True):
        self.subject_path = subject_path
        db_types = {
            "nucl": ["blastn", "tblastn"],
            "prot": ["blastp", "blastx"]
        }
        self.db_type = db_type.lower()
        if self.db_type not in db_types:
            exit("Error: Unrecognised BLAST database type:", self.db_type)
        index_suffices_suffices = ["hr", "in", "sq", "db", "ot", "tf", "to"]
        self.index_suffices = {
            "nucl": [".n" + x for x in index_suffices_suffices],
            "prot": [".p" + x for x in index_suffices_suffices]
        }
        self.program_type = program_type.lower()
        if program_type not in db_types[self.db_type]:
            exit("Error: Unrecognised BLAST program name:", program_type, "for database type:", self.db_type)
        if do_indexing:
            self.makeIndex()
        self.searchprocess = None
        self.max_hits = max_hits
        # search queue and results for batch searches
        self.batch_search_queue = []
        self.batch_search_result = {}
    
    def makeIndex(self):
        """
        Index subject_path using makeblastdb.
        """
        index_found = True
        for suffix in self.index_suffices[self.db_type]:
            index_found = index_found and os.path.exists(self.subject_path + suffix)
        if not index_found:
            subprocess.run(["makeblastdb", "-dbtype", self.db_type, "-in", self.subject_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE, text=True)
    
    def reIndex(self):
        """
        Remove and remake index files for subject_path.
        """
        self.removeIndex()
        self.makeIndex()

    def removeIndex(self):
        """
        Remove index files for subject_path.
        """
        for suffix in self.index_suffices[self.db_type]:
            os.remove(self.subject_path + suffix)
    
    def sortBlastHits(self, hits, sortby="bitscore", reverse=True):
        """
        Sort BLAST hits by a given attribute.
        """
        if sortby not in ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]:
            exit("Error: Unrecognised BLAST hit attribute:", sortby)
        return sorted(hits, key=lambda x: getattr(x, sortby), reverse=reverse)
    
    def doSearch(self, queries):
        """
        Do a search using a list of queries (with .fasta() and .sequenceName() methods) against subject_path.
        """
        self.searchprocess = subprocess.Popen([self.program_type, "-db", self.subject_path, "-outfmt", "6", "-num_threads", str(multiprocessing.cpu_count())], stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE, text=True)
        stdout_lines = self.searchprocess.communicate(input="\n".join([query.fasta() for query in queries]))[0].splitlines()
        # limit hits if max_hits set
        if self.max_hits is not None:
            stdout_lines = stdout_lines[:self.max_hits]
        return self.sortBlastHits([BlastHit(x) for x in stdout_lines])
    
    def queueBatchSearch(self, queries):
        """
        Queue a search using queries against subject_path.
        """
        # queue queries if not already queued or in results
        for query in queries:
            if query not in self.batch_search_queue and query.sequenceName() not in self.batch_search_result:
                self.batch_search_queue.append(query)
    
    def getBatchSearch(self, queries):
        """
        Get the results of a batch search, returned as a dict indexed by query sequence name.
        """
        # queue queries, whcih checks for duplication
        self.queueBatchSearch(queries)
        # trigger search if there are unprocessed queued queries
        if len(self.batch_search_queue) > 0:
            self.triggerBatchSearch()
        # get results and return as dict by sequence name
        query_names = [query.sequenceName() for query in queries]
        results = {}
        for query_name in query_names:
            if query_name in self.batch_search_result:
                results[query_name] = self.batch_search_result[query_name]
            else:
                results[query_name] = []
        return results

    def triggerBatchSearch(self):
        """
        Trigger a batch search of the queued queries and append to batch_search_result.
        """
        # do search and reset queue
        results = self.doSearch(self.batch_search_queue)
        self.batch_search_queue = []
        # append results to batch_search_result
        for result in results:
            if result.qseqid not in self.batch_search_result:
                self.batch_search_result[result.qseqid] = []
            self.batch_search_result[result.qseqid].append(result)
    
    def clearBatchSearch(self):
        """
        Clear the batch search queue and result.
        """
        self.batch_search_queue = []
        self.batch_search_result = {}

class Orf:
    """
    ORF object, used by NuclSequence, representing a single open reading frame in a DNA sequence with a maximal extent start and stop and alternative start codons.
    Works with the forward strand, so reverse complement the sequence if necessary and coordinates are in the reverse direction.
    """
    def __init__(self, sequence: Union["NuclSequence", "ProtSequence"], start:int, stop:int, forward:bool, altstarts:list=[]):
        self.sequence = sequence # parent chromosome/contig/mRNA sequence
        self.start = start
        self.stop = stop
        self.forward = forward # True if forward strand, False if reverse strand
        self.length = self.stop - self.start
        self.altstarts = sorted(altstarts)
        self.frame = stop % 3
    
    @property
    def name(self):
        return self.sequence.sequence_name + "_" + str(self.start) + "-" + str(self.stop) + ("f" if self.forward else "r")

    def __repr__(self):
        return fasta(self)
    
    def dnaseq(self):
        return NuclSequence(self.sequence.forward[self.start:self.stop])

    def protseq(self):
        return self.dnaseq().protein().sequence
    
    def sequenceName(self):
        return self.sequence.sequence_name + "_" + str(self.start) + "-" + str(self.stop) + "-" + ("f" if self.forward else "r")

    def fasta(self):
        return ">" + self.sequenceName() + "\n" + self.protseq()
    
    def orfLengthPvalue(self):
        """
        Simple p value of an ORF of up to or equal to that length (nucleotides), assuming random 50% GC sequence
        """
        return (1 - len(self.sequence.codons["*"]) / self.dnaseq().ncodons) ** (self.length / 3)
    
    def orfStartPvalue(self):
        """
        Simple p value of an ORF having started by that distance through the sequence, assuming random 50% GC sequence
        """
        return (1 - len(self.sequence.codons["M"]) / self.dnaseq().ncodons) ** self.start

    def findOrfsFromStop(self):
        """
        Working from the stop codon, longest possible orf and alternative start codons
        """
        offset = 3
        starts = []
        # walk left from the stop codon, until another stop codon found or the end of the sequence reached
        while self.stop-offset-3 >=0 and self.sequence.forward[self.stop-offset-3:self.stop-offset] not in self.sequence.codons["*"] and "N" not in self.sequence.forward[self.stop-offset-3:self.stop-offset]:
            # record start codons found on the way
            if self.sequence.forward[self.stop-offset-3:self.stop-offset] in self.sequence.codons["M"]:
                starts.append(offset)
            offset += 3
        if len(starts) > 0:
            # take the last found start codon as start, also records alternative starts
            # [doing necessary coordinate adjustment]
            self.start = self.stop - starts[-1] - 3
            self.altstarts = sorted([self.stop - x - 3 for x in starts])
        else:
            # if no start codons found, set start to stop, no alternative starts
            self.start = self.stop
            self.altstarts = []
        self.length = self.stop - self.start

    def selectLongest(self):
        """
        Set the start codon to the longest ORF defined by altstarts
        """
        self.start = self.altstarts[-1]
        self.length = self.stop - self.start

    def blastTestStart(self, blast):
        """
        Test all start codons in altstarts for a dropoff in tBLASTn bitscore, and set the start cdon to the first non-outlier
        Score each start-start and final start-stop codon segment by bitscore per amino acid
        Extend rightwards to minquerylength amino acids for short segments
        """
        def stdev(list):
            list = [x for x in list if x is not None]
            return statistics.stdev(list) if len(list) > 1 else 0
        def mean(list):
            list = [x for x in list if x is not None]
            return statistics.mean(list) if len(list) > 0 else None

        thresholdstdevs = 1.5 # threshold for outlier detection
        minquerylength = 25 * 3 # 25 amino acids, 3 nucleotides per codon
        zerohitsbitscore = 0 # bitscore to impute for no BLAST hits
        # make query segments starting at each start codon
        rawranges = [x - self.start for x in self.altstarts + [self.stop]]
        # derive start/end ranges, extend rightward to minimum query length if too short
        ranges = []
        for index in range(len(rawranges) - 1):
            start = rawranges[index]
            end = rawranges[index + 1]
            if rawranges[index + 1] - rawranges[index] < minquerylength:
                end = min(rawranges[index] + minquerylength, rawranges[-1])
            ranges.append((start, end))
        scores = []
        # do blast search for each segment
        seqs = []
        for index in range(len(ranges)):
            seqs.append(ProtSequence(self.dnaseq().protein().sequence[ranges[index][0]//3:ranges[index][1]//3], sequence_name=self.sequenceName() + "_startindex" + str(index)))
        for index in range(len(ranges)):
            result = blast.getBatchSearch([seqs[index]])[seqs[index].sequenceName()]
            score = zerohitsbitscore
            if len(result) > 0:
                score = result[0].bitscore / seqs[index].length() # normalise by length
            scores.append(score)
        # search for first non-outlier start codon
        mean, stdev = mean(scores), stdev(scores)
        threhsold = mean - thresholdstdevs * stdev
        index = 0
        while scores[index] < threhsold and index < len(scores) - 1:
            index += 1
        changed = self.altstarts[index] != self.start
        # set output
        self.start = self.altstarts[index] if index < len(self.altstarts) else self.stop
        self.length = self.stop - self.start
        # return report
        #if changed: print([round(x, 3) for x in [mean, stdev, threhsold]], [round(x, 3) for x in scores])
        return changed, index

class NuclSequence:
    """
    DNA sequence object for handling a single DNA sequence, including a naive three-frame ORF finder
    """
    def __init__(self, sequence:str, sequence_name:str="unnamed"):
        self.forward = sequence.upper()
        self.sequence_name = sequence_name
        self.alphabet={"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
        self.codons = {
            "A": ["GCT", "GCC", "GCA", "GCG"],
            "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
            "N": ["AAT", "AAC"],
            "D": ["GAT", "GAC"],
            "C": ["TGT", "TGC"],
            "Q": ["CAA", "CAG"],
            "E": ["GAA", "GAG"],
            "G": ["GGT", "GGC", "GGA", "GGG"],
            "H": ["CAT", "CAC"],
            "I": ["ATT", "ATC", "ATA"],
            "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
            "K": ["AAA", "AAG"],
            "F": ["TTT", "TTC"],
            "P": ["CCT", "CCC", "CCA", "CCG"],
            "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
            "T": ["ACT", "ACC", "ACA", "ACG"],
            "Y": ["TAT", "TAC"],
            "V": ["GTT", "GTC", "GTA", "GTG"],
            "M": ["ATG"],
            "W": ["TGG"],
            "*": ["TAA", "TAG", "TGA"]
        }
        self.ncodons = 0
        self.aas = {}
        for aa in self.codons:
            self.ncodons += len(self.codons[aa])
            for v in self.codons[aa]:
                self.aas[v] = aa
    
    def __repr__(self):
        return self.fasta()
    
    def __len__(self):
        return len(self.forward)
    
    def length(self):
        return self.__len__()

    def sequenceName(self):
        return self.sequence_name

    def fasta(self, case:str="upper"):
        if case == "upper":
            return ">" + self.sequence_name + "\n" + self.forward
        else:
            return ">" + self.sequence_name + "\n" + self.forward.lower()

    def _reverseComplement(self, sequence):
        return "".join([self.alphabet[x] for x in sequence[::-1]])

    def reverseComplement(self):
        return self._reverseComplement(self.forward)

    def _translation(self, sequence):
        if "N" not in sequence:
            return "".join([self.aas[x] for x in [sequence[i:i+3] for i in range(0, len(sequence)-len(sequence)%3, 3)]])
        return None

    def protein(self):
        return ProtSequence(self._translation(self.forward), sequence_name=self.sequence_name)

    def startCodon(self):
        return self.forward[:3] in self.codons["M"]
    
    def stopCodon(self):
        return self.forward[-3:] in self.codons["*"]
    
    def containsNs(self):
        return "N" in self.forward
    
    def inFrame(self):
        return len(self.forward) % 3 == 0
    
    def isOrf(self):
        return self.startCodon() and self.stopCodon() and not self.containsNs() and self.inFrame()
    
    def findOrfs(self, min_length:int=150, search_rev_compl:bool=False):
        """
        Finds ORFs by finding stop codons and finding the furthest upstream start codon without an intervening in-frame stop
        Only returns ORFs over min_length bases (min_length / 3 amino acids) long
        """
        def searchStrand(sequence, min_length, is_rev_compl=False):
            # search reverse complement for reverse strand ORFs, resulting coordinates are on reverse sequence
            if is_rev_compl:
                forward = sequence.reverseComplement()
            else:
                forward = sequence.forward
            # find stops
            stops = []
            for i in range(len(forward) - 3):
                if forward[i:i+3] in self.codons["*"]:
                    stops.append(i+3)
            # find longest orf for each stop, and all possible alternate starts
            hits = []
            for stop in stops:
                # for each stop, set up a stub of an ORF from the stop codon
                hit = Orf(NuclSequence(forward, sequence_name=sequence.sequence_name), stop, stop, not is_rev_compl, altstarts=[])
                # find all possible start codons upstream from the stop
                hit.findOrfsFromStop()
                if hit.length > min_length:
                    # if sufficiently long (no upstream start codon gives length zero), add to hits
                    hits.append(hit)
            return hits
        
        hits = searchStrand(self, min_length=min_length, is_rev_compl=False)
        if search_rev_compl:
            hits += searchStrand(self, min_length=min_length, is_rev_compl=True)
        return hits

class ProtSequence:
    """
    Protein sequence object for handling a single protein sequence.
    """
    def __init__(self, sequence:str, sequence_name:str="unnamed"):
        self.sequence = sequence.upper()
        self.sequence_name = sequence_name
        self.alphabet = "ACDEFGHIKLMNPQRSTVWY*X"

    def __repr__(self):
        return self.fasta()

    def __len__(self):
        return len(self.sequence)

    def length(self):
        return self.__len__()

    def sequenceName(self):
        return self.sequence_name

    def fasta(self, case:str="upper"):
        if case == "upper":
            return ">" + self.sequence_name + "\n" + self.sequence
        else:
            return ">" + self.sequence_name + "\n" + self.sequence.lower()

class OrfFinder:
    """
    For running OrfFinder on a query FASTA file
    """
    def __init__(self, reference_genomes_path:str, query_fasta_path:str, verbosity:int=2):
        self.query_fasta = Fasta(path=query_fasta_path)
        self.blast = Blast(reference_genomes_path, "nucl", "tblastn")
        self.verbosity = verbosity
    
    def allGoodOrfs(self, min_length:int=300, threshold_bit_score:float=35, search_rev_compl:bool=True):
        """
        Returns the best ORFs in the query fasta, preserving long ORFs, removing overlapping and low tBLASTn bitscore ORFs.
        Checks the remaining ORFs for the best start codon using tBLASTn.
        Suited for finding all ORFs in a genome.
        """
        result = {}
        for sequence in self.query_fasta.sequences:
            if self.verbosity > 1: print("Sequence:", sequence)
            orfs = self.query_fasta.sequences[sequence].findOrfs(min_length=min_length, search_rev_compl=search_rev_compl)
            if self.verbosity > 1: print("", "Open reading frames:", len(orfs), "over", min_length, "bp")
            # first, remove overlaps, iterating from the longest ORFs, removing shorter ones in any overlap
            orfs = sorted(orfs, key=lambda x: x.length, reverse=True)
            filteredorfs = orfs.copy()
            for i in range(len(orfs)):
                left, right = min(orfs[i].start, orfs[i].stop), max(orfs[i].start, orfs[i].stop)
                for j in range(i+1, len(orfs)):
                    if left < orfs[j].stop < right or left < orfs[j].start < right:
                        filteredorfs[j] = None
            orfs = [x for x in filteredorfs if x is not None]
            if self.verbosity > 1: print("", "Longest non-overlapping:", len(orfs))
            # second, remove ORFs with low blast bitscore (< threshold_bit_score)
            filteredorfs = [x for x in orfs if x.protseq() is not None]
            for orf in orfs:
                # queue query for batch search
                self.blast.queueBatchSearch([orf])
            for orf in orfs:
                # get batch search results and parse
                blastresult = self.blast.getBatchSearch([orf])[orf.sequenceName()]
                bitscore = blastresult[0].bitscore if len(blastresult) > 0 else 0
                if bitscore < threshold_bit_score:
                    filteredorfs.remove(orf)
            orfs = filteredorfs
            if self.verbosity > 1: print("", "Passed tBLASTn filtering:", len(orfs))
            # finally, adjust start codon for each ORF
            if self.verbosity > 1: print("", "Adjusting start codon by tBLASTn checks:")
            for orf in orfs:
                changed, index = orf.blastTestStart(self.blast)
                if self.verbosity > 1 and changed: print("", "", "Changed start codon:", orf.sequenceName(), "to methionine", index + 1)
            if len(orfs) > 0:
                result[sequence] = orfs
        return result
    
    def bestOrf(self, min_length:int=300):
        """
        Returns the longest, leftmost (first) and best tBLASTn hit ORF for each sequence in the query FASTA.
        Checks the best hit for the best start codon using tBLASTn.
        Suited for finding the best ORF in a transcript.
        """
        def orfByName(orfs:list, name:str):
            for orf in orfs:
                if orf.sequenceName() == name:
                    return orf
            return None

        result = {}
        for sequence in self.query_fasta.sequences:
            if self.verbosity > 1: print("Sequence:", sequence)
            orfs = self.query_fasta.sequences[sequence].findOrfs(min_length=min_length, search_rev_compl=False)
            if self.verbosity > 1: print("", "Open reading frames:", len(orfs), "over", min_length, "bp")
            if len(orfs) > 0:
                longest = sorted(orfs, key=lambda x: x.length, reverse=True)[0] # longest ORF found
                first = sorted(orfs, key=lambda x: x.start, reverse=False)[0] # leftmost ORF found
                blasthits = self.blast.doSearch("\n".join(x.fasta() for x in orfs)) # tBLASTn all orfs to find best hit
                blastbitscore = blasthits[0].bitscore if len(blasthits) > 0 else 0
                blastevalue = blasthits[0].evalue if len(blasthits) > 0 else None
                blasthit = orfByName(orfs, blasthits[0].qseqid) if len(blasthits) > 0 else None # top bitscore-sorted tBLASTn hit
                if self.verbosity > 1:
                    print("", "Best ORFs for", sequence)
                    print("", "", "Longest ORF:", longest.sequenceName(), longest.length, "bp long", "e=", longest.orfLengthPvalue())
                    print("", "", "First ORF:", first.sequenceName(), first.start, "bp from sequence start", "e=", first.orfStartPvalue())
                    print("", "", "Best tBLASTn hit:", blasthit.sequenceName(), blastbitscore, "bitscore", "e=", blastevalue)
                hits = list(dict.fromkeys([x for x in [longest, first, blasthit] if x is not None]))
                # adjust start codon for each unique hit
                if self.verbosity > 1: print("", "Adjusting start codon by tBLASTn checks:")
                for hit in hits:
                    changed, index = hit.blastTestStart(self.blast)
                    if self.verbosity > 1 and changed: print("", "", "Changed start codon:", hit.sequenceName(), "to methionine", index + 1)
                result[sequence] = hits
            else:
                if self.verbosity > 1: print("No ORFs found for", sequence)
        return result
    
    def adjustOrfs(self):
        """
        Adjusts start codons for all ORFs in the query FASTA based on tBLASTn hits.
        """
        result = {}
        for sequence in self.query_fasta.sequences:
            # TODO, handle loading of existing ORFs from GFF file
            # orfs = [loading function]
            exit("Not implemented!!!")
            for orf in orfs:
                changed, index = orf.blastTestStart(self.blast)
                if self.verbosity > 1 and changed: print("", "", "Changed start codon:", orf.sequenceName(), "to methionine index", index)
            if len(orfs) > 0:
                result[sequence] = orfs
        return result

class Fasta:
    """
    Fasta object, handling FASTA file reading and writing, or interpreting a string as FASTA.
    .sequences is a dict of NuclSequence objects, indexed by name.
    """
    def __init__(self, path:str=None, string:str=None):
        # string, or read string from file
        self.string = None
        if string is not None:
            self.string = string
        elif path is not None:
            with open(path, "r") as file:
                self.string = file.read()
        # sequences as a dict of NuclSequence objects, indexed by name
        self.sequences = None
        if self.string is not None:
            self.sequences = {x.splitlines()[0].split(" ")[0]: NuclSequence("".join(x.splitlines()[1:]).upper(), sequence_name=x.splitlines()[0].split(" ")[0]) for x in self.string.split(">")[1:]}

    def __repr__(self):
        return fasta(self)

    def fasta(self):
        return "\n".join([">" + x[0] + "\n" + x[1] for x in self.sequences.items()])
