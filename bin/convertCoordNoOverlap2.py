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
  parser.add_argument('-f','--maxOverlap',dest='max_acceptable_overlap',default=0.1,help='Max percentage of sequence overlap')
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
              domainevalue=vec[11] #e-value for domain
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

def load_annotation_pfam2(filename,reference_dict):
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
              domainevalue=vec[11] #e-value for domain
              try: #make sure coordinates can be converted to int
                  startcoord=int(vec[19]) #internal coordinates
                  stopcoord=int(vec[20]) #internal coordinates
              except ValueError:
                     print "coordinates or PFAM file in wrong format"
                     raise
              #convert coordinates
              (newcontig,newstart,newstop,strandinfo)=convert_coordinates_one_seq(reference_dict,contigID,startcoord,stopcoord)
              
              #check if domain exists in dictionary and add fields
              if newcontig not in annot.keys():
                 annot[newcontig]=[]
                 start[newcontig]=[]
                 stop[newcontig]=[]
                 strand[newcontig]=[]
                 evalue[newcontig]=[]
                 annot[newcontig].append(domain)
                 start[newcontig].append(newstart)
                 stop[newcontig].append(newstop)
                 strand[newcontig].append(strandinfo)
                 evalue[newcontig].append(domainevalue)
              else:
                    annot[newcontig].append(domain)
                    start[newcontig].append(newstart)
                    stop[newcontig].append(newstop)
                    strand[newcontig].append(strandinfo)
                    evalue[newcontig].append(domainevalue)
    f.close()
    print "Number of contigs loaded: "+str(len(annot.keys()))
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
              line=newcontig+'\t'+str(newstart)+'\t'+str(newstop)+'\t'+strand+'\t'+domain+'\n'
            else:
              line=newcontig+'\tTIGRFAM\tprotein_hmm_match\t'+str(newstart)+'\t'+str(newstop)+'\t.\t'+strand+'\t.\tID='+domain+'\n'
            g.write(line)
    g.close()
 
def convert_coordinates_one_seq(reference_dict,contig,start,stop):
    '''function for converting coordinates from amino acid to nucleotide
    and saving in Tentacle output format'''

    record_dict=reference_dict
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
    return (newcontig,newstart,newstop,strand)
    
def report_only_best_overlapping_sequence(annot,startvec,stopvec,strand,evalue,maxAcceptableOverlap):
  newannot={}
  newstartvec={}
  newstopvec={}
  newstrand={}
  newevalue={}
  for contig in annot:
    besthits={}
    besthits['+']={}
    besthits['-']={}
    
    if len(annot[contig])>1:
      for i in range(len(annot[contig])):
        strandi=strand[contig][i]
        newhit=True
        for hit in besthits[strandi]:
          startj=int(hit.split('-')[0])
          stopj=int(hit.split('-')[1])
          overlap=min(stopvec[contig][i],stopj)-max(startvec[contig][i],startj)+1
          minlen=min(stopvec[contig][i]-startvec[contig][i]+1,stopj-startj+1)
          maxoverlap=minlen*maxAcceptableOverlap
          if overlap > maxoverlap: #overlapping, save the best one
            if float(evalue[contig][i])<float(evalue[contig][besthits[strandi][hit]]): #best hit
              newkey=str(startvec[contig][i])+'-'+str(stopvec[contig][i])
              besthits[strandi][newkey]=i
              del besthits[strandi][hit]
              newhit=False
              break #save each hit max once
            else:
              newhit=False
        if newhit==True:
          newkey=str(startvec[contig][i])+'-'+str(stopvec[contig][i])
          besthits[strandi][newkey]=i
    else: #only one hit
      newkey=str(startvec[contig][0])+'-'+str(stopvec[contig][0])
      besthits[strand[contig][0]][newkey]=0
    newannot[contig]=[]
    newstartvec[contig]=[]
    newstopvec[contig]=[]
    newstrand[contig]=[]
    newevalue[contig]=[]
    for strandi in besthits:
      for hit in besthits[strandi]:
        hitindex=besthits[strandi][hit]
        newannot[contig].append(annot[contig][hitindex])
        newstartvec[contig].append(startvec[contig][hitindex])
        newstopvec[contig].append(stopvec[contig][hitindex])
        newstrand[contig].append(strand[contig][hitindex])
        newevalue[contig].append(evalue[contig][hitindex])
  return (newannot,newstartvec,newstopvec,newstrand,newevalue)

def writeToFile(output_filename,tentacle_format,annot,startvec,stopvec,strand,evalue):
  g=open(output_filename,'w')
  if not tentacle_format:
    g.write('##gff-version 3\n')
  for contig in annot:
    for i in range(len(annot[contig])):
      if tentacle_format:
        #tentacle file format
        line=contig+'\t'+str(startvec[contig][i])+'\t'+str(stopvec[contig][i])+'\t'+strand[contig][i]+'\t'+annot[contig][i]+'\t'+evalue[contig][i]+'\n'
      else:
        line=contig+'\tTIGRFAM\tprotein_hmm_match\t'+str(startvec[contig][i])+'\t'+str(stopvec[contig][i])+'\t.\t'+strand[contig][i]+'\t.\tID='+annot[contig][i]+'\n'
      g.write(line)
  g.close()

    
def main(reference_filename,annot_filename,output_filename,tentacle_format,max_acceptable_overlap):
    record_dict=index_fasta(reference_filename)
    (annot,startvec,stopvec,strand,evalue)=load_annotation_pfam2(annot_filename,record_dict)
    (annot,startvec,stopvec,strand,evalue)=report_only_best_overlapping_sequence(annot,startvec,stopvec,strand,evalue,max_acceptable_overlap)
    writeToFile(output_filename,tentacle_format,annot,startvec,stopvec,strand,evalue)
         
    
    
    
    
    
    
    
if __name__=='__main__':
  arguments=parseArgsCC()
  main(arguments.reference_file,arguments.annotation_file,arguments.output_file,arguments.tentacle_format,arguments.max_acceptable_overlap)