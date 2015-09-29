#!/usr/bin/python

def parseArgs(argv):
''' Function for parsing arguments 
'''

  import argparse
  parser = argparse.ArgumentParser(description='Extract PFAM or TIGRFAM domains from one or several fasta files and create one fasta file per domain as output')
  parser.add_argument('-r', dest='reference_dir',help='The path to the reference directory')
  parser.add_argument('-a', dest='annotation_dir', help='The path to the reference directory')
  parser.add_argument('-o', dest='output_dir', help='The path to the output directory')
  parser.add_argument('-n', dest='n',help='number of threads')
  arguments=parser.parse_args(sys.argv[1:])
  
def main(r,a,o,n):
  
  

if __name__=='__main__':
  import sys
  arguments=parseArgs(sys.argv[1])
  main(arguments.reference_dir,arguments.annotation_dir,arguments.output_dir,arguments.n)