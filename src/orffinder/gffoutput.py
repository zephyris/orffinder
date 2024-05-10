import geffa

from .orffinder import Orf

class GffOutput:
    def __init__(self, query_fasta_file):
        self.gff: geffa.GffFile = geffa.GffFile(fasta_file=query_fasta_file)
    
    def add_orf(self, orf: Orf):
        sequence_region = self.gff.sequence_regions[orf.sequence.sequencename]
        if orf.forward:
            gff_start = orf.start + 1
            gff_end = orf.stop + 1
            gff_strand = '+'
        else:
            sequence_length = orf.sequence.length()
            gff_start = sequence_length - orf.stop
            gff_end = sequence_length - orf.start
            gff_strand = '-'
        gene_name = orf.orfname()
        geffa.geffa.GeneNode(-1, sequence_region, "orffinder", 'gene', gff_start, gff_end, '0.0', gff_strand, '.', f'ID={gene_name}')
        geffa.geffa.MRNANode(-1, sequence_region, "orffinder", 'mRNA', gff_start, gff_end, '0.0', gff_strand, '.', f'ID={gene_name}.1;Parent={gene_name}')
        geffa.geffa.CDSNode(-1, sequence_region, "orffinder", 'CDS', gff_start, gff_end, '0.0', gff_strand, '0', f'ID={gene_name}.1-CDS1;Parent={gene_name}.1')

    def save(self):
        return str(self.gff)