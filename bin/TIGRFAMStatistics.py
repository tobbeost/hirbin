#!/usr/bin/python
def getTIGRFAMstruct(path):
    import glob
    filelist=glob.glob(path+'*.hmmout')
    annot={}
    for filename in filelist:
        #  function for loading annotation and coordinates from COG, Pfam and TIGRfam
        import os
        #list the annotation files in annotation directory
        print "loading annotation from file " + filename
        
    
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


                  
                  #check if domain exists in dictionary and add fields
                  if domain not in annot.keys():
                     annot[domain]=[]
               
                     annot[domain].append(contigID)
 
                  else:
                        annot[domain].append(contigID)
        f.close()
        print "Number of motifs loaded: "+str(len(annot.keys()))
        #return annotation information as a tuple.
    output=annot
    return(output)   
    

    
def main(path):
    annot=getTIGRFAMstruct(path)
    g=open('numberofcontigsTIGRFAM.txt','w')
    for domain in annot:
        g.write(domain+'\t'+str(len(annot[domain]))+'\n')
    g.close() 
    
    
    
    
    
    
    

    
if __name__=="__main__":
   import sys
   path=sys.argv[1]
   main(path)