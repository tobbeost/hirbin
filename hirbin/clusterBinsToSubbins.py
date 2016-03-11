#!/usr/bin/env python2.6
# coding: utf-8
#HirBin function clusterBinsToSubbins
from parsers import *
import os
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
import argparse
import sys

def parseArgs():
  ''' 
    Function for parsing arguments 
  '''
  parser = argparse.ArgumentParser(description='Convert coordinates from protein sequence to nucleotide sequence coordinates.')
  parser.add_argument('-m', dest='mapping_file',help='The path to the mapping/metadata file (required)',required=True)
  parser.add_argument('-o', dest='output_dir', help='The name of the output directory. If not specified, a new directory is created')
  parser.add_argument('-id' '--sequenceIdentity', dest='identity',default=0.7,help='Sequence identity cutoff used (default:  %(default)s)')
  arguments=parser.parse_args(sys.argv[1:]) 
  return arguments

def createOutputDirectory(output_directory):
  ''' Function for creating a new output directory. If the name is not specified, decide a new name '''
  if output_directory==None: #no output directory specified, create a new output directory
    new_output_dir='hirbin_output'
    not_created_dir=False
    suffix=1
    #create a new name for the output directory
    while not not_created_dir:
      if os.path.isdir(new_output_dir):
        suffix=suffix+1
        new_output_dir='hirbin_output_'+str(suffix)
      else:
        not_created_dir=True
    try:
      os.mkdir(new_output_dir)
      output_directory=new_output_dir
      print "Creating a new output directory at "+os.getcwd()+'/'+output_directory
    except OSError as e:
      if not force:
        raise
  else:
    if not os.path.isdir(output_directory):
      try:
        os.mkdir(output_directory)
      except OSError as e:
        raise
  print output_directory
  return(output_directory)

def main(mappingFile,output_directory,identity):
  
  output_directory=createOutputDirectory(output_directory)
  print "Reading metadata file"
  metadata=Hirbin_run(output_directory)
  metadata.readMetadata(mappingFile)
  print metadata.groups
  print metadata.reference
  print identity

if __name__=='__main__':
  arguments=parseArgs()
  main(arguments.mapping_file,arguments.output_dir,arguments.identity)
