import os, sys, subprocess, statistics, multiprocessing

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
        # search queue for batch searches
        self.batch_search_queue = []
        self.batch_search_result = {}
    
    def makeIndex(self):
        """
        Index subject_path using makeblastdb
        """
        index_found = True
        for suffix in self.index_suffices[self.db_type]:
            index_found = index_found and os.path.exists(self.subject_path + suffix)
        if not index_found:
            subprocess.run(["makeblastdb", "-dbtype", self.db_type, "-in", self.subject_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE, text=True)
    
    def reIndex(self):
        """
        Remove and remake index files for subject_path
        """
        self.removeIndex()
        self.makeIndex()

    def removeIndex(self):
        """
        Remove index files for subject_path
        """
        for suffix in self.index_suffices[self.db_type]:
            os.remove(self.subject_path + suffix)
    
    def sortBlastHits(self, hits, sortby="bitscore", reverse=True):
        """
        Sort BLAST hits by a given attribute
        """
        if sortby not in ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]:
            exit("Error: Unrecognised BLAST hit attribute:", sortby)
        return sorted(hits, key=lambda x: getattr(x, sortby), reverse=reverse)
    
    def doSearch(self, query_lines):
        """
        Do a search using query against subject_path
        """
        self.searchprocess = subprocess.Popen([self.program_type, "-db", self.subject_path, "-outfmt", "6", "-num_threads", str(multiprocessing.cpu_count())], stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE, text=True)
        stdout_lines = self.searchprocess.communicate(input=query_lines)[0].splitlines()
        if self.max_hits is not None:
            stdout_lines = stdout_lines[:self.max_hits]
        return self.sortBlastHits([BlastHit(x) for x in stdout_lines])
    
    def queueBatchSearch(self, query_lines):
        """
        Queue a search using query against subject_path, indexed in queue by query id
        """
        query_fasta = Fasta(string=query_lines)
        for query_id in query_fasta.sequences:
            if query_id not in self.batch_search_queue and query_id not in self.batch_search_result:
                self.batch_search_queue.append(query_fasta.sequences[query_id].fasta())
    
    def getBatchSearch(self, query_lines):
        """
        Get the results of a batch search, returned as a dict indexed by query id
        """
        self.queueBatchSearch(query_lines)
        if len(self.batch_search_queue) > 0:
            self.triggerBatchSearch()
        query_ids = Fasta(string=query_lines).sequences.keys()
        results = {}
        for query_id in query_ids:
            if query_id in self.batch_search_result:
                results[query_id] = self.batch_search_result[query_id]
            else:
                results[query_id] = []
        return results

    def triggerBatchSearch(self):
        """
        Trigger a batch search of the queued queries and append to batch_search_result
        """
        results = self.doSearch("\n".join(self.batch_search_queue)+"\n")
        self.batch_search_queue = []
        for result in results:
            if result.qseqid not in self.batch_search_result:
                self.batch_search_result[result.qseqid] = []
            self.batch_search_result[result.qseqid].append(result)
    
    def clearBatchSearch(self):
        """
        Clear the batch search queue and result
        """
        self.batch_search_queue = []
        self.batch_search_result = {}

class Orf:
    """
    ORF object, used by DnaSequence, representing a single open reading frame in a DNA sequence with a maximal extent start and stop and alternative start codons
    Works with the forward strand, so reverse complement the sequence if necessary
    """
    def __init__(self, sequence, start, stop, forward, altstarts=[]):
        self.sequence = sequence # parent chromosome/contig/mRNA sequence
        self.start = start
        self.stop = stop
        self.forward = forward # True if forward strand, False if reverse strand
        self.length = self.stop - self.start
        self.altstarts = sorted(altstarts)
        self.frame = stop % 3
    
    def __repr__(self):
        return "\t".join(str(x) for x in [self.orfname(), self.protseq()])
    
    def dnaseq(self):
        return DnaSequence(self.sequence.forward[self.start:self.stop])

    def protseq(self):
        return self.dnaseq().translation()
    
    def orfname(self):
        return self.sequence.sequencename + "_" + str(self.start) + "-" + str(self.stop) + "-" + ("f" if self.forward else "r")

    def fasta(self):
        return "\n".join([">" + self.orfname(), self.protseq(), ""])
    
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
        for index in range(len(ranges)):#
            seqs.append(DnaSequence(self.dnaseq().forward[ranges[index][0]:ranges[index][1]], sequencename=self.orfname() + "_startindex" + str(index)))
            blast.queueBatchSearch(seqs[-1].fastaProt())
        for index in range(len(ranges)):
            result = blast.getBatchSearch(seqs[index].fastaProt())[seqs[index].sequencename]
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

class DnaSequence:
    """
    DNA sequence object, including a naive three-frame ORF finder
    """
    def __init__(self, sequence, sequencename="unnamed"):
        self.forward = sequence.upper()
        self.sequencename = sequencename
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
        return self.forward
    
    def _reverseComplement(self, sequence):
        return "".join([self.alphabet[x] for x in sequence[::-1]])
    
    def reverseComplement(self):
        return self._reverseComplement(self.forward)
    
    def translation(self):
        if "N" not in self.forward:
            return "".join([self.aas[x] for x in [self.forward[i:i+3] for i in range(0, len(self.forward)-len(self.forward)%3, 3)]])
        return None
    
    def fasta(self):
        return ">" + self.sequencename + "\n" + self.forward
    
    def fastaProt(self):
        return ">" + self.sequencename + "\n" + self.translation()
    
    def startCodon(self):
        return self.forward[:3] in self.codons["M"]
    
    def stopCodon(self):
        return self.forward[-3:] in self.codons["*"]
    
    def containsNs(self):
        return "N" in self.forward
    
    def length(self):
        return len(self.forward)

    def __len__(self):
        return self.length()
    
    def findOrfs(self, minlength = 150, searchrevcompl = False):
        """
        Finds ORFs by finding stop codons and finding the furthest upstream start codon without an intervening in-frame stop
        Only returns ORFs over minlength bases (minlength / 3 amino acids) long
        """        
        def searchStrand(sequence, minlength, revcompl=False):
            # search reverse complement for reverse strand ORFs, resulting coordinates are on reverse sequence
            if revcompl:
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
                hit = Orf(DnaSequence(forward, sequencename=sequence.sequencename), stop, stop, not revcompl, altstarts=[])
                # find all possible start codons upstream from the stop
                hit.findOrfsFromStop()
                if hit.length > minlength:
                    # if sufficiently long (no upstream start codon gives length zero), add to hits
                    hits.append(hit)
            return hits
        
        hits = searchStrand(self, minlength, revcompl=False)
        if searchrevcompl:
            hits += searchStrand(self, minlength, revcompl=True)
        return hits

class OrfFinder:
    def __init__(self, reference_genomes_path, query_fasta_path, verbosity=2):
        self.query_fasta = Fasta(path=query_fasta_path)
        self.blast = Blast(reference_genomes_path, "nucl", "tblastn")
        self.verbosity = verbosity
    
    def allGoodOrfs(self, minlength=300, thresholdbitscore=35, searchrevcompl=True):
        """
        Returns the best ORFs in the query fasta, preserving long ORFs, removing overlapping and low tBLASTn bitscore ORFs.
        Checks the remaining ORFs for the best start codon using tBLASTn.
        Suited for finding all ORFs in a genome.
        """
        result = {}
        for sequence in self.query_fasta.sequences:
            if self.verbosity > 1: print("Sequence:", sequence)
            orfs = self.query_fasta.sequences[sequence].findOrfs(minlength=minlength, searchrevcompl=searchrevcompl)
            if self.verbosity > 1: print("", "Open reading frames:", len(orfs), "over", minlength, "bp")
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
            # second, remove ORFs with low blast bitscore (< thresholdbitscore)
            filteredorfs = [x for x in orfs if x.protseq() is not None]
            for orf in orfs:
                # queue query for batch search
                self.blast.queueBatchSearch(orf.fasta())
            for orf in orfs:
                # get batch search results and parse
                blastresult = self.blast.getBatchSearch(orf.fasta())[orf.orfname()]
                bitscore = blastresult[0].bitscore if len(blastresult) > 0 else 0
                if bitscore < thresholdbitscore:
                    filteredorfs.remove(orf)
            orfs = filteredorfs
            if self.verbosity > 1: print("", "Passed tBLASTn filtering:", len(orfs))
            # finally, adjust start codon for each ORF
            if self.verbosity > 1: print("", "Adjusting start codon by tBLASTn checks:")
            for orf in orfs:
                changed, index = orf.blastTestStart(self.blast)
                if self.verbosity > 1 and changed: print("", "", "Changed start codon:", orf.orfname(), "to methionine", index + 1)
            if len(orfs) > 0:
                result[sequence] = orfs
        return result
    
    def bestOrf(self, minlength=300):
        """
        Returns the longest, leftmost (first) and best tBLASTn hit ORF for each sequence in the query FASTA.
        Checks the best hit for the best start codon using tBLASTn.
        Suited for finding the best ORF in a transcript.
        """
        def orfByName(orfs, name):
            for orf in orfs:
                if orf.orfname() == name:
                    return orf
            return None

        result = {}
        for sequence in self.query_fasta.sequences:
            if self.verbosity > 1: print("Sequence:", sequence)
            orfs = self.query_fasta.sequences[sequence].findOrfs(minlength=minlength, searchrevcompl=False)
            if self.verbosity > 1: print("", "Open reading frames:", len(orfs), "over", minlength, "bp")
            if len(orfs) > 0:
                longest = sorted(orfs, key=lambda x: x.length, reverse=True)[0] # longest ORF found
                first = sorted(orfs, key=lambda x: x.start, reverse=False)[0] # leftmost ORF found
                blasthits = self.blast.doSearch("\n".join(x.fasta() for x in orfs)) # tBLASTn all orfs to find best hit
                blastbitscore = blasthits[0].bitscore if len(blasthits) > 0 else 0
                blastevalue = blasthits[0].evalue if len(blasthits) > 0 else None
                blasthit = orfByName(orfs, blasthits[0].qseqid) if len(blasthits) > 0 else None # top bitscore-sorted tBLASTn hit
                if self.verbosity > 1:
                    print("", "Best ORFs for", sequence)
                    print("", "", "Longest ORF:", longest.orfname(), longest.length, "bp long", "e=", longest.orfLengthPvalue())
                    print("", "", "First ORF:", first.orfname(), first.start, "bp from sequence start", "e=", first.orfStartPvalue())
                    print("", "", "Best tBLASTn hit:", blasthit.orfname(), blastbitscore, "bitscore", "e=", blastevalue)
                hits = list(dict.fromkeys([x for x in [longest, first, blasthit] if x is not None]))
                # adjust start codon for each unique hit
                if self.verbosity > 1: print("", "Adjusting start codon by tBLASTn checks:")
                for hit in hits:
                    changed, index = hit.blastTestStart(self.blast)
                    if self.verbosity > 1 and changed: print("", "", "Changed start codon:", hit.orfname(), "to methionine", index + 1)
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
                if self.verbosity > 1 and changed: print("", "", "Changed start codon:", orf.orfname(), "to methionine index", index)
            if len(orfs) > 0:
                result[sequence] = orfs
        return result

class Fasta:
    def __init__(self, path=None, string=None):
        self.string = None
        if string is not None:
            self.string = string
        elif path is not None:
            with open(path, "r") as file:
                self.string = file.read()
        self.sequences = None
        if self.string is not None:
            self.sequences = {x.splitlines()[0].split(" ")[0]: DnaSequence("".join(x.splitlines()[1:]), sequencename=x.splitlines()[0].split(" ")[0]) for x in self.string.split(">")[1:]}
    
    def fasta(self):
        return "\n".join([">" + x[0] + "\n" + x[1] for x in self.sequences.items()])

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
