#!/usr/bin/env python
# coding: utf-8
import os
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
import argparse
import sys

def parseArgsCC():
  ''' 
    Function for parsing arguments 
  '''
  parser = argparse.ArgumentParser(description='Convert coordinates from protein sequence to nucleotide sequence coordinates.')
  parser.add_argument('-r', dest='reference_file',help='The path to the assembly contigs/reference file (required)',required=True)
  parser.add_argument('-a', dest='annotation_file', help='The path to the annotation file (required)',required=True)
  parser.add_argument('-o', dest='output_file', help='The name of the output file',required=True)
  parser.add_argument('--tentacle_format', dest='tentacle_format', action='store_true',help='Write the output in Tentacle format (default: gff format)')
  arguments=parser.parse_args(sys.argv[1:]) 
  return arguments

def load_annotation_pfam(filename):
    '''function for loading annotation and coordinates from PFAM nd TIGRfam
    (HMMersearch domtblout format) and saving into dictionaries'''

    #list the annotation files in annotation directory
    print "loading annotation from file " + filename
    annot={}
    start={}
    stop={}
    strand={}
    evalue={}
    #open the file
    f=open(filename)
    #skip header lines
    f.readline()
    f.readline()
    f.readline()
    
    #extract the domain and sequence_id information from annotation file
    for line in f:
        line=line.rstrip('\n')
        vec=line.split()
        if len(vec)>20: #extra check
           domain=vec[3] #domain ID
           if domain not in "NA":
              contigID=vec[0] #sequence id
              domainevalue=vec[11] #e-value
              try: #make sure coordinates can be converted to int
                  startcoord=int(vec[19]) #internal coordinates
                  stopcoord=int(vec[20]) #internal coordinates
              except ValueError:
                     print "coordinates or PFAM file in wrong format"
                     raise
              
              strandinfo='+'    #all features are on sense strand
              #check if domain exists in dictionary and add fields
              if domain not in annot.keys():
                 annot[domain]=[]
                 start[domain]=[]
                 stop[domain]=[]
                 strand[domain]=[]
                 evalue[domain]=[]
                 annot[domain].append(contigID)
                 start[domain].append(startcoord)
                 stop[domain].append(stopcoord)
                 strand[domain].append(strandinfo)
                 evalue[domain].append(domainevalue)
              else:
                    annot[domain].append(contigID)
                    start[domain].append(startcoord)
                    stop[domain].append(stopcoord)
                    strand[domain].append(strandinfo)
                    evalue[domain].append(domainevalue)
    f.close()
    print "Number of motifs loaded: "+str(len(annot.keys()))
    #return annotation information as a tuple.
    output= (annot,start,stop,strand,evalue)
    return(output)   


def index_fasta(fasta_path):
    '''function for indexing the amino acid fasta reference file'''
    print "No fasta index file found, indexing reference fasta file"
    
    record_dict = SeqIO.index(fasta_path, "fasta")
    return(record_dict)

    
def convert_coordinates(reference_dict,annot,startvec,stopvec,strand,evalue,output_filename,tentacle_format):
    '''function for converting coordinates from amino acid to nucleotide
    and saving in Tentacle output format'''

    record_dict=reference_dict
    g=open(output_filename,'w')
    if not tentacle_format:
      g.write('##gff-version 3\n')
    for domain in annot.keys():
        
        for i in range(len(annot[domain])):
            contig=annot[domain][i]
            try:
                start=startvec[domain][i]
                stop=stopvec[domain][i]
            except ValueError:
                   print "wrong format"
                   raise
            contig2=contig.split('_')
            frame=contig2.pop()
            newcontig='_'.join(contig2)
           
            if frame in '1':
               newstart=start*3
               newstop=stop*3
               strand='+'
            elif frame in '2':
                 newstart=start*3+1
                 newstop=stop*3+1
                 strand='+'
            elif frame in '3':
                 newstart=start*3+2
                 newstop=stop*3+2
                 strand='+'
            elif frame in '4':
                 n=len(str(record_dict[contig].seq))
                 newstart=n*3-(stop*3)+1
                 newstop=n*3-(start*3)+1
                 strand='-'
            elif frame in '5':
                 n=len(str(record_dict[contig].seq))
                 newstart=n*3-(stop*3)  +2
                 newstop=n*3-(start*3)+2
                 strand='-'
            elif frame in '6':
                 n=len(str(record_dict[contig].seq))
                 newstart=n*3-(stop*3)+3
                 newstop=n*3-(start*3)+3
                 strand='-'
            if tentacle_format:
              #tentacle file format
              line=newcontig+'\t'+str(newstart)+'\t'+str(newstop)+'\t'+strand+'\t'+domain+'\t'+evalue[domain][i]+'\n'
            else:
              line=newcontig+'\tTIGRFAM\tprotein_hmm_match\t'+str(newstart)+'\t'+str(newstop)+'\t.\t'+strand+'\t.\tID='+domain+'\n'
            g.write(line)
    g.close()
    
    
def main(reference_filename,annot_filename,output_filename,tentacle_format):
    record_dict=index_fasta(reference_filename)
    (annot,startvec,stopvec,strand,evalue)=load_annotation_pfam(annot_filename)
    convert_coordinates(record_dict,annot,startvec,stopvec,strand,evalue,output_filename,tentacle_format)
         
    
    
    
    
    
    
    
if __name__=='__main__':
  arguments=parseArgsCC()
  main(arguments.reference_file,arguments.annotation_file,arguments.output_file,arguments.tentacle_format)