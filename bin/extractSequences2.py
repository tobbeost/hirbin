#!/usr/bin/env python
#one sample at a time for parallelization purposes
def load_annotation_pfam(filename):
    #  function for loading annotation and coordinates from COG, Pfam and TIGRfam
    #input file= .hmm file from hmmsearch
    import os
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
    #samplename=filename.rstrip('.seq.fa.pep.hmmout')
    #extract the domain and sequence_id information from annotation file
    for line in f:
        line=line.rstrip('\n')
        vec=line.split()
        if len(vec)>20: #extra check
           domain=vec[3] #domain ID
           if domain not in "NA":
              contigID=vec[0] #sequence id
              domainevalue=vec[6] #e-value
              try: #make sure coordinates can be converted to int
                  startcoord=int(vec[19]) #internal coordinates
                  stopcoord=int(vec[20]) #internal coordinates
              except ValueError:
                     print "coordinates or PFAM file in wrong format"
                     break
              
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
    print "indexing fasta file " + fasta_path
    from Bio import SeqIO
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    from Bio.SeqRecord import SeqRecord
    record_dict = SeqIO.index(fasta_path, "fasta")
    return(record_dict)

def extract_sequences(args):
    (fasta_path,annotation_path, output_dir)=args
    domainInfo=load_annotation_pfam(annotation_path)
    sample=fasta_path.rstrip('.fasta').split('/')[-1]
    print "Generating domain fasta sequences for "+sample+" ..."
    from Bio import SeqIO
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    from Bio.SeqRecord import SeqRecord 
    (annot,start,stop,strand,evalue)=domainInfo
    record_dict=index_fasta(fasta_path)
    recordlist=[]
    outfilename=output_dir + sample +'.fasta'
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


def main(reference_dir,annotation_path,output_dir,n):
    import os, os.path
    try:
        n=int(n)
    except ValueError:
           print "n should be integer"
    
    from multiprocessing.dummy import Pool
    fasta_list=os.listdir(reference_dir)
    annot_list=os.listdir(annotation_path)
    fasta_list.sort()
    annot_list.sort()
    arglist=[(reference_dir+i,annotation_path+j,output_dir) for (i,j) in zip(fasta_list,annot_list)]
    print "Initiating " + str(n) + " processes"
    p=Pool(n)
    p.map(extract_sequences,arglist)
    
   
if __name__=='__main__':
   import argparse, os, sys
   parser = argparse.ArgumentParser(description='Extract PFAM or TIGRFAM domains from one or several fasta files and create one fasta file per domain as output')
   parser.add_argument('-r', dest='reference_dir',help='The path to the reference directory')
   parser.add_argument('-a', dest='annotation_dir', help='The path to the reference directory')
   parser.add_argument('-o', dest='output_dir', help='The path to the output directory')
   parser.add_argument('-n', dest='n',help='number of threads')
   arguments=parser.parse_args(sys.argv[1:])
   main(arguments.reference_dir,arguments.annotation_dir,arguments.output_dir,arguments.n)