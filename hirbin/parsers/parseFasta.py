#!/usr/bin/env python
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
def index_fasta(fasta_path):
    print "indexing fasta file " + fasta_path
    record_dict = SeqIO.index(fasta_path, "fasta")
    return(record_dict)