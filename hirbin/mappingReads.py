#!/usr/bin/env python
from hirbin import *
from hirbin.parsers import *
#from parsers import *
import os
from os import mkdir
import argparse
import sys
from multiprocessing.dummy import Pool

def parseArgs(argv):
  ''' 
    Function for parsing arguments 
  '''
  parser = argparse.ArgumentParser(description='Perform mapping of reads to annotated reference sequences, by running bowtie2 with default parameters.')
  parser.add_argument('-m', dest='mapping_file',help='The path to the mapping/metadata file (required)',required=True)
  parser.add_argument('-o', dest='output_dir', help='The path to the output directory. If not specified, a new directory is created')
  parser.add_argument('-n', dest='n',help='number of threads [default = %(default)s]',default=1, type=int)
  #parser.add_argument('-f',dest='force_output_directory',action="store_true",help='Force, if the output directory should be overwritten')
  arguments=parser.parse_args(argv)
  return arguments


def indexReferenceFile(referenceFile,sampleName,outputdir):
  command="bowtie2-build "+referenceFile+" "+sampleName+" > "+ outputdir+"/bowtieoutput/bowtie_log_"+sampleName+".txt 2> "+outputdir+"/bowtieoutput/bowtie_stats_"+sampleName+".txt"
  os.system(command)

def mapping(fastq1,fastq2,sampleName,outputdir):
  command="bowtie2 --no-unal -x "+sampleName+" -1 "+fastq1+" -2 "+fastq2+" -S "+outputdir+"/bowtieoutput/"+sampleName+".sam >> "+ outputdir+"/bowtieoutput/bowtie_log_"+sampleName+".txt 2>> "+outputdir+"/bowtieoutput/bowtie_stats_"+sampleName+".txt"
  os.system(command)
  
def convertToBam(referenceFile,sampleName,outputdir):
  command="samtools import "+referenceFile+" "+outputdir+"/bowtieoutput/"+sampleName+".sam "+outputdir+"/bowtieoutput/"+sampleName+".bam"
  os.system(command)
  os.system("rm "+outputdir+"/bowtieoutput/"+sampleName+".sam")

def calculateCoverage(sampleName,outputdir):
  command="bedtools coverage -abam "+outputdir+"/bowtieoutput/"+sampleName+".bam -b "+outputdir+"/"+sampleName+".gff -counts > "+outputdir+"/bowtieoutput/"+sampleName+".cov"
  os.system(command)
  parseCoverageBed(outputdir+"/bowtieoutput/"+sampleName+".cov",outputdir+"/bowtieoutput/"+sampleName+".tab")
  
def runMapping(metadata,ncpus):
  p=Pool(ncpus)
  output_dir=metadata.output_directory
  arglist=zip(metadata.reference.keys(),metadata.reference.values(),metadata.reads1.values(),metadata.reads2.values(),[output_dir]*len(metadata.reference.keys()))
  p.map(mappingOneSample,arglist)
  
  
def mappingOneSample(args):
  (sampleName,referenceFile,fastq1,fastq2,outputdir)=args
  print "Running Bowtie2 for sample "+sampleName
  indexReferenceFile(referenceFile,sampleName,outputdir)
  mapping(fastq1,fastq2,sampleName,outputdir)
  convertToBam(referenceFile,sampleName,outputdir)
  calculateCoverage(sampleName,outputdir)
  
def main(mapping_file,output_directory,ncpus,force):
  metadata=Hirbin_run(output_directory)
  metadata.readMetadata(mapping_file)
  if not force:
    output_directory=metadata.createOutputDirectory(output_directory)
  metadata.output_directory=output_directory
  try:
    mkdir(output_directory+'/bowtieoutput/')
  except OSError as e:
    if not force:
      print "Output directory already exists, you can use an already existing output directory by including the flag -f"
      raise
  runMapping(metadata,ncpus)
  
if __name__=='__main__':
  args=parseArgs(sys.argv[1:])
  force=True
  main(args.mapping_file,args.output_dir,args.n,force)

