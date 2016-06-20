from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

def load_annotation_pfam(filename):
    '''function for loading annotation and coordinates from PFAM and TIGRfam
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
              domainevalue=vec[6] #e-value
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