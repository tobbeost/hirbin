#!/usr/bin/env python
# coding: utf-8

import hirbin.parsers
import os
from os import path, mkdir
from os.path import isdir
import glob
import argparse
import sys
from multiprocessing.dummy import Pool


def parseArgs(argv):
  ''' 
    Function for parsing arguments 
  '''
  parser = argparse.ArgumentParser(description='Perform functional annotation of assembled metagenomics sequences, using TIGRFAM or PFAM.')
  parser.add_argument('-m', dest='mapping_file',help='The path to the mapping/metadata file (required)',required=True)
  parser.add_argument('-db', dest='database_dir', help='The path to the protein family database directory (required)',required=True)
  parser.add_argument('--type', dest='type', default='nucleotide', help='Sequence type. "prot" if the reference fasta files are for amino acid sequences. "nucl" if nucleotide sequences. In the latter case the nucleotide sequences will be translated. (default:  %(default)s)')
  parser.add_argument('-o', dest='output_dir', help='The path to the output directory. If not specified, a new directory is created')
  parser.add_argument('-e','--e-value', dest='evalue_cutoff', default='1e-10', help='E-value cutoff [default = %(default)s]')
  parser.add_argument('-n', dest='n',help='number of threads [default = %(default)s]',default=1, type=int)
  parser.add_argument('-f',dest='force_output_directory',action="store_true",help='Force, if the output directory should be overwritten')
  arguments=parser.parse_args(argv)
  return arguments


  
def runTranslation(sample_dict,ncpus,output_dir):
  ''' Translate the nucleotide sequences in the assembly w.r.t all 6 reading frames using transeq, using n cpu:s'''
  p=Pool(ncpus)
  arglist=[path+' '+output_dir+'/protseq/'+name+'.pep' for (name,path) in sample_dict.iteritems()]
  p.map(translateOneSample,arglist)
  return(output_dir+'/protseq/')
  
def translateOneSample(argument):
  '''Subroutine for translating one sample on one cpu'''
  os.system('transeq --frame 6 ' + argument)
  
def runHMMer(sample_dict,ncpus,output_dir,database_dir,evalue_cutoff,protseq_dir):
  ''' Run HMMer to perform functional annotation for all samples using n cpu:s'''
  p=Pool(ncpus)
  arglist=['--domtblout ' + output_dir +'/hmmeroutput/' + name + '.hmmout -E ' + evalue_cutoff + ' ' + database_dir + ' ' +  protseq_dir + '/'+ name +'.pep' for name in sample_dict.keys()]
  print arglist
  p.map(functionalAnnotationOneSample,arglist)

def functionalAnnotationOneSample(argument):
  ''' Run HMMer for one sample using 1 cpu'''
  os.system('hmmsearch '+ argument +'> /dev/null') #dispose unwanted output to /dev/null

def main(mapping_file,database_dir,output_directory,type,ncpus,evalue_cutoff,force):
  metadata=Hirbin_run(output_directory)
  metadata.readMetadata(mapping_file)
  output_directory=metadata.createOutputDirectory(output_directory)
  metadata.output_directory=output_directory
  if type.startswith("nucl"):
    #if protein sequences are missing, run translation
    try:
      mkdir(output_directory+'/protseq/')
    except OSError as e:
      if not force:
        raise
    protseq_dir=runTranslation(metadata.reference,ncpus,output_directory)
    try:
      mkdir(output_directory+'/hmmeroutput/')
    except OSError as e:
      if not force:
        raise
    runHMMer(metadata.reference,ncpus,output_directory,database_dir,evalue_cutoff,protseq_dir)
    #perform functional annotation using HMMer

if __name__=='__main__':
  arguments=parseArgs(sys.argv[1:])
  main(arguments.mapping_file,arguments.database_dir,arguments.output_dir,arguments.type,arguments.n,arguments.evalue_cutoff,arguments.force_output_directory)
