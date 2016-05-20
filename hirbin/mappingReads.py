#!/usr/bin/env python
import hirbin.parsers
import os
import argparse
import sys

def parseArgs(argv):
  ''' 
    Function for parsing arguments 
  '''
  parser = argparse.ArgumentParser(description='Perform mapping of reads to annotated reference sequences.')
  parser.add_argument('-m', dest='mapping_file',help='The path to the mapping/metadata file (required)',required=True)
  parser.add_argument('-db', dest='database_dir', help='The path to the protein family database directory (required)',required=True)
  parser.add_argument('--type', dest='type', default='nucleotide', help='Sequence type. "prot" if the reference fasta files are for amino acid sequences. "nucl" if nucleotide sequences. In the latter case the nucleotide sequences will be translated. (default:  %(default)s)')
  parser.add_argument('-o', dest='output_dir', help='The path to the output directory. If not specified, a new directory is created')
  parser.add_argument('-e','--e-value', dest='evalue_cutoff', default='1e-10', help='E-value cutoff [default = %(default)s]')
  parser.add_argument('-n', dest='n',help='number of threads [default = %(default)s]',default=1, type=int)
  parser.add_argument('-f',dest='force_output_directory',action="store_true",help='Force, if the output directory should be overwritten')
  arguments=parser.parse_args(argv)
  return arguments


  
def mapping()

