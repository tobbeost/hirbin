def load_annotation_pfam(filename):
    '''function for loading annotation and coordinates from PFAM nd TIGRfam
    (HMMersearch domtblout format) and saving into dictionaries'''
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


def index_fasta(fasta_path,sample):
    '''function for indexing the amino acid fasta reference file'''
    print "No fasta index file found, indexing reference fasta file"
    import marshal
    from Bio import SeqIO
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    from Bio.SeqRecord import SeqRecord
    record_dict = SeqIO.index(fasta_path, "fasta")
    return(record_dict)
    c=open(sample+'.fasta.index','wb')
    marshal.dump(record_dict,c)
    c.close()
    
def convert_coordinates(reference_dict,annot,startvec,stopvec,strand,evalue,output_filename):
    '''function for converting coordinates from amino acid to nucleotide
    and saving in Tentacle output format'''

    record_dict=reference_dict
    g=open(output_filename,'w')
    for domain in annot.keys():
        
        for i in range(len(annot[domain])):
            contig=annot[domain][i]
            try:
                start=startvec[domain][i]
                stop=stopvec[domain][i]
            except ValueError:
                   print "wrong format"
                   break
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
            line=newcontig+'\t'+str(newstart)+'\t'+str(newstop)+'\t'+strand+'\t'+domain+'\n'
            g.write(line)
    g.close()
    
    
def main(sample):
    
    reference_filename='/storage/tobiaso/metagenomes/qin2012/contigs/protseq/'+sample+'.contigs.m100.gz.pep'
    annot_filename='TIGRFAM/'+sample+'.contigs.m100.gz.hmmout'
    output_filename='TIGRFAM3/'+sample+'.annotation.tab'
    import cPickle,os
    from Bio import SeqIO
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    from Bio.SeqRecord import SeqRecord
    path='TIGRFAM/'
    #if os.path.isfile(sample+'.fasta.index'):
    #   c=open(sample+'.fasta.index','rb')
    #   record_dict=cPickle.load(c)
    #else:
    record_dict=index_fasta(reference_filename,sample)
    #print "converting annotation"
    
    #filename=sample+'.fa.hmmout'
    (annot,startvec,stopvec,strand,evalue)=load_annotation_pfam(annot_filename)  
    convert_coordinates(record_dict,annot,startvec,stopvec,strand,evalue,output_filename)
         
    
    
    
    
    
    
    
if __name__=='__main__':
   import sys
   sample=sys.argv[1]
   main(sample)