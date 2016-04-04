#!/usr/bin/env python
# coding: utf-8
#HirBin function clusterBinsToSubbins

from hirbin/parsers import *
from hirbin/thirdparty import runUclust
import os
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
import argparse
import sys
import re
from multiprocessing.dummy import Pool

def parseArgs():
  ''' 
    Function for parsing arguments 
  '''
  parser = argparse.ArgumentParser(description='Cluster annotated sequences (bins) into sub-bins')
  parser.add_argument('-m', dest='mapping_file',help='The path to the mapping/metadata file (required)',required=True)
  parser.add_argument('--type', dest='type', default='protein', help='Sequence type. "prot" if the reference and annotation are for amino acid sequences. "nucl" if nucleotide sequences. In the latter case the nucleotide sequences will be translated. (default:  %(default)s)')
  parser.add_argument('-o', dest='output_dir', help='The name of the output directory. If not specified, a new directory is created')
  parser.add_argument('-id', '--sequenceIdentity', dest='identity',default=0.7,help='Sequence identity cutoff used for sub-binning (default:  %(default)s)')
  parser.add_argument('-p','--minRepresented', dest='p', default=0.75, help='Require non-zero counts in at least the fraction p of the samples for the bin/sub-bin to be representative. (default:  %(default)s)')
  parser.add_argument('--minMeanCount', dest='minMeanCount', default=3, help='Minimum mean count per sample for the bin/sub-bin to be representative. (default:  %(default)s)')
  parser.add_argument('-n', dest='n',help='number of threads (default = %(default)s)',default=1, type=int)
  parser.add_argument('-f',dest='force_output_directory',action="store_true",help='Force, if the output directory should be overwritten')
  arguments=parser.parse_args(sys.argv[1:]) 
  return arguments

  
def extract_sequences_one_sample(args):
    (fasta_path,annotation_path,sample, output_dir)=args
    domainInfo=load_annotation_pfam(annotation_path)
    print "Generating domain fasta sequences for "+sample+" ..."
    from Bio import SeqIO
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    from Bio.SeqRecord import SeqRecord 
    (annot,start,stop,strand,evalue)=domainInfo
    record_dict=index_fasta(fasta_path)
    recordlist=[]
    outfilename=output_dir +'/forClustering/'+ sample +'.fasta'
    outhandle=open(outfilename,'w')
    for domainID in annot.keys():
        for i in range(len(annot[domainID])):
            domain=annot[domainID][i]
            try:
                seq=record_dict[domain]
            except KeyError:
               print "Error: " + domain + " not in fasta file.\n"
               break    
            a=start[domainID][i]
            b=stop[domainID][i]
            seq_strand=strand[domainID][i]
            seq_evalue=evalue[domainID][i]
            if seq_strand in '+':
               domain_feature = SeqFeature(FeatureLocation(a-1, b-1), type="domain", strand=1)
            elif seq_strand in '-':
                 domain_feature = SeqFeature(FeatureLocation(a-1, b-1), type="domain", strand=-1)
            feature_seq = domain_feature.extract(seq)
            feature_seq.id=feature_seq.id+' '+domainID+' '+seq_evalue 
            recordlist.append(feature_seq)
    
    
    SeqIO.write(recordlist, outhandle, "fasta")
    outhandle.close()
    
    print "done"

def getSequencesPerDomain(metadata):
  outputdir=metadata.output_directory
  filelist=[outputdir+'/forClustering/'+samplename+".fasta" for samplename in metadata.reference.keys()]
  domains={}
  for filename in filelist:
    f1=open(filename)
    sequences = SeqIO.parse(f1, "fasta")
    for s in sequences:
      domain=s.description.split()[1]
      if domain not in domains.keys():
        domains[domain]=[]
        domains[domain].append(s)
      else:
        domains[domain].append(s)
    f1.close()
  
  for domain in domains:
    g=open(outputdir+'/forClustering/'+domain+'.fasta','w')
    SeqIO.write(domains[domain], g, "fasta")
    g.close()
  for filename in filelist:
    os.remove(filename)


def extractSequences(metadata,n,force):
  outputdir=metadata.output_directory
  try:
    os.mkdir(outputdir+'/forClustering/')
  except OSError as e:
    if not force:
      print "Directory "+outputdir+'/forClustering/ already exists'
      raise
  fasta_list=metadata.reference.values()
  annot_list=metadata.annotation.values()
  sample_list=metadata.annotation.keys()
  arglist=[(i,j,sample,outputdir) for (i,j,sample) in zip(fasta_list,annot_list,sample_list)]
  print "Initiating " + str(n) + " processes"
  p=Pool(n)
  p.map(extract_sequences_one_sample,arglist)
  getSequencesPerDomain(metadata)
  
def getSubBins(groups,clustpath,cutoffnumber,minMeanCount,identity,countDict):
  samplelisttotal=groups.keys()
  cutoff=round(len(samplelisttotal)*cutoffnumber)
  domainlist=os.listdir(clustpath)
  directory="clust"+str(identity) #loop through each ID cutoff
  print "getting subBins for " + directory
  countsvec={}#initiate counts vector
  countsvec[directory]={} #initiate counts vector
  for fname in domainlist:
    domain=fname.rstrip('.fasta.uc')
    countsvec[directory][domain]={} #initiate counts vector
    clusters=getClusterStruct(clustpath+'/'+fname)
    for key in clusters.keys(): #for each cluster, check if the cluster is large enough for statistical testing
      #if n>=cutoff: #if the cluster is large enough, get the counts from results file.
      samplelist=[]
      n=0
      for value in clusters[key]:
        contigID=re.sub('_[0-9]$','',value) #remove trailing number
        try:
          crow=countDict[contigID+'_'+domain]
      
          for sample in crow:
            if sample not in samplelist and countDict[contigID+'_'+domain][sample]>0:
              samplelist.append(sample)
              n+=1
        except KeyError:
          print "No count data fouond for"+ contigID+'_'+domain
          raise
             
      if n>=cutoff:
        countsvec[directory][domain][key]={}#initiate counts vector
        totcount=0
        for value in clusters[key]:
          contigID=re.sub('_[0-9]$','',value) #remove trailing number
          try:
            crow=countDict[contigID+'_'+domain]
            for sample in crow:
              if sample in countsvec[directory][domain][key].keys():
                countsvec[directory][domain][key][sample]+=countDict[contigID+'_'+domain][sample]
              else:
                countsvec[directory][domain][key][sample]=countDict[contigID+'_'+domain][sample]
                  #check mean count
          except KeyError:
            print "No count data fouond for"+ contigID+'_'+domain
            raise
        meancount=sum(countsvec[directory][domain][key].values())/float(len(samplelist))
        if meancount<minMeanCount: #remove if mean count is less than minMeanCount.
          del(countsvec[directory][domain][key])
             
  
  #countsvec[directory][domain][key][sample]=samplecount
  
  print "writing counts to file for " + directory
  g=open(clustpath+'/../abundance_matrix_subbins_'+directory+'.txt','w')
  samplelisttotal=sorted(samplelisttotal)
  for s in samplelisttotal:
      g.write('\t'+s)
  g.write('\n')
  for domain in countsvec[directory].keys():
      for key in countsvec[directory][domain].keys():
          g.write(domain+'_'+directory+'_'+str(key))
          for s in samplelisttotal:
              if s in countsvec[directory][domain][key]:
                 g.write('\t'+str(countsvec[directory][domain][key][s]))
              else:
                   g.write('\t0')
          g.write('\n')
  g.close()
  print "done"



def main(mappingFile,output_directory,type,p,minMeanCount,identity,n,force):
  #reading metadata file and creating output directory
  metadata=Hirbin_run(output_directory)
  metadata.readMetadata(mappingFile)
  output_directory=metadata.createOutputDirectory(output_directory)
  metadata.output_directory=output_directory
  extractSequences(metadata,n,force)
  try:
    os.mkdir(output_directory+'/clust'+str(identity))
  except OSError as e:
    if not force:
      raise
  runUclust(output_directory+"/forClustering/",identity)
  countDict=getCountStruct(metadata)
  getSubBins(metadata.groups,output_directory+"/clust"+str(identity),p,minMeanCount,identity,countDict)
  domains=createAbundanceMatrix(metadata,p,minMeanCount)

if __name__=='__main__':
  arguments=parseArgs()
  main(arguments.mapping_file,arguments.output_dir,arguments.type,arguments.p,arguments.minMeanCount,arguments.identity,arguments.n,arguments.force_output_directory)
