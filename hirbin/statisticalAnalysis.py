#!/usr/bin/env python
# coding: utf-8
#HirBin function statisticalAnalysis
from hirbin import *
from hirbin.parsers import *
import os
import argparse
import sys
import pkg_resources
def parseArgs():
  ''' 
    Function for parsing arguments 
  '''
  parser = argparse.ArgumentParser(description='Run statistical analysis')
  parser.add_argument('-m', dest='mapping_file',help='The path to the mapping/metadata file (required)',required=True)
  parser.add_argument('-o', dest='output_dir', help='The name of the hierbin output directory, created in a previous run.',required=True)
  parser.add_argument('-ref', '--referenceGroup', dest='ref',help='Reference/control condition used for calculating log fold change. If not specified, the last group listed in the metadata file is used as reference/control.')
  arguments=parser.parse_args(sys.argv[1:]) 
  return arguments

def main(mapping_file,output_dir,ref):
  metadata=Hirbin_run(output_dir)
  metadata.readMetadata(mapping_file)
  output_directory=metadata.createOutputDirectory(output_dir)
  metadata.output_directory=output_dir
  groups=[]
  for sample in metadata.samples:
    groups.append(metadata.groups[sample])
  groupslist=','.join(groups)
  if (ref is None):
    ref=groups[-1]
  filelist=os.listdir(output_dir)
  for filename in filelist:
    if filename.startswith('abundance_matrix'):
      filename=output_dir+'/'+filename
      print "Running statistical analysis for "+filename
      command=groupslist+' '+ref+' '+filename
      path = pkg_resources.resource_filename('hirbin', 'statistical_analysis.R')
      os.system('Rscript '+path+ ' ' + command)
  
if __name__=='__main__':
  arguments=parseArgs()
  main(arguments.mapping_file,arguments.output_dir,arguments.ref)
