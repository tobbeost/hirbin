#!/usr/bin/env python
# coding: utf-8
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
  parser.add_argument('-o', dest='output_dir', help='The path to the output directory')
  parser.add_argument('-e','--e-value', dest='evalue_cutoff', default='1e-10', help='E-value cutoff [default = %(default)s]')
  parser.add_argument('-p',dest='protein_seq_dir',help='directory of protein sequences. If omitted it will be created by translating the nucleotide sequences.')
  parser.add_argument('-n', dest='n',help='number of threads [default = %(default)s]',default=1, type=int)
  parser.add_argument('-f',dest='force_output_directory',action="store_true",help='Force, if the output directory should be overwritten')
  arguments=parser.parse_args(argv)
  #if arguments.output_dir==None: #no output directory specified, create a new output directory
  #if arguments.mapping_file==None:
  
  return arguments


def main(mapping_file,database_dir,output_directory,ncpus,evalue_cutoff,force,protseq_dir):
  output_directory=createOutputDirectory(output_directory,protseq_dir,force)
  sample_dict=readMappingFile(mapping_file)
  if protseq_dir==None: #if protein sequences are missing, run translation
    protseq_dir=runTranslation(sample_dict,ncpus,output_directory) 
  runHMMer(sample_dict,ncpus,output_directory,database_dir,evalue_cutoff,protseq_dir) #perform functional annotation using HMMer
  
def runTranslation(sample_dict,ncpus,output_dir):
  ''' Translate the nucleotide sequences in the assembly w.r.t all 6 reading frames using transeq, using n cpu:s'''
  p=Pool(ncpus)
  arglist=[path+' '+output_dir+'/protseq/'+name+'.pep' for (name,path) in sample_dict.iteritems()]
  p.map(translateOneSample,arglist)
  return('protseq/')
  
def translateOneSample(argument):
  '''Subroutine for translating one sample on one cpu'''
  os.system('transeq --frame 6 ' + argument)
  
def runHMMer(sample_dict,ncpus,output_dir,database_dir,evalue_cutoff,protseq_dir):
  ''' Run HMMer to perform functional annotation for all samples using n cpu:s'''
  p=Pool(ncpus)
  arglist=['--domtblout ' + output_dir +'/' + name + '.hmmout -E ' + evalue_cutoff + ' ' + database_dir + ' ' + output_dir +'/'+ protseq_dir + name +'.pep' for name in sample_dict.keys()]
  print arglist
  p.map(functionalAnnotationOneSample,arglist)

def functionalAnnotationOneSample(argument):
  ''' Run HMMer for one sample using 1 cpu'''
  os.system('hmmsearch '+ argument +'> /dev/null') #dispose unwanted output to /dev/null
  

def readMappingFile(mapping_file):
  '''
  Reads in the sample mapping/metadata file and creates a dictionary 
  with sample names and file path to assembly for that sample.
  '''
  with open(mapping_file) as f:
    line=f.readline() #skip header
    samples={}
    for line in f:
      line=line.rstrip()
      line=line.split('\t')
      samplename=line[0]
      samplepath=line[-1]
      if samplename not in samples:
        samples[samplename]=samplepath
    return samples

def createOutputDirectory(output_directory,protseq_dir,force):
  ''' Function for creating a new output directory. If the name is not specified, decide a new name '''
  if output_directory==None: #no output directory specified, create a new output directory
    new_output_dir='hierbin_output'
    not_created_dir=False
    suffix=1
    #create a new name for the output directory
    while not not_created_dir:
      if os.path.isdir(new_output_dir):
        suffix=suffix+1
        new_output_dir='hierbin_output_'+str(suffix)
      else:
        not_created_dir=True
    try:
      mkdir(new_output_dir)
      if protseq_dir==None:
        mkdir(new_output_dir+'/protseq/')
      output_directory=new_output_dir
    except OSError as e:
      if not force:
        raise
  else:
    try:
      mkdir(output_directory)
      if protseq_dir==None:
        mkdir(output_directory+'/protseq/')
    except OSError as e:
      if not force:
        raise
  print output_directory
  return(output_directory)
                
if __name__=='__main__':
  arguments=parseArgs(sys.argv[1:])
  main(arguments.mapping_file,arguments.database_dir,arguments.output_dir,arguments.n,arguments.evalue_cutoff,arguments.force_output_directory,arguments.protein_seq_dir)
