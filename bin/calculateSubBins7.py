#!/usr/bin/python
def getClusterStruct(filename):
    #function for reading one cluster structure for one PFAM domain and save it into dictonary
    #filename, file name of cluster file
    #returns a dictonary with contigs/genes belonging to each cluster.
    domainID=filename.split('/')[-1].split('.')[0]
    f=open(filename)
    clusters={}
    sum=0
    for line in f:
        if line[0] in ["S","H"] :
           line=line.split("\t")
           clusterID=int(line[1])
           target=line[8]
           if clusterID not in clusters:
              clusters[clusterID]=[]
              clusters[clusterID].append(target)
           else:
                clusters[clusterID].append(target)
       
    f.close()
    return clusters
    
def getTIGRFAMstruct(filelist):
    for filename in filelist:
        #  function for loading annotation and coordinates from COG, Pfam and TIGRfam
        import os
        #list the annotation files in annotation directory
        print "loading annotation from file " + filename
        annot={}
    
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


def getCountStruct(results):
    #function for reading all results from mapping and save in a dictonary
    print "Reading and saving mapping results..."
    import glob, re
    filelist=glob.glob(results+'*.fastq.tab')
    countDict={}
    n=len(filelist)
    r=0
    q=0
    test=range(n/10,n,n/10)
    for filename in filelist:
        sample=filename.split('/')[-1].split('.')[0]
        r=r+1
        if r in test:
           q+=1
           print repr(q*10) + '%... '
        f=open(filename)
        for line in f:
            parts=line.split('\t')
            contigID=parts[0].split(':')[0]
            counts=int(parts[1])
            if contigID not in countDict:
               countDict[contigID]={}
            if sample not in countDict[contigID]:
               countDict[contigID][sample]=counts
            else:
                 countDict[contigID][sample]+=counts
        f.close()
    return countDict

def getSubBins(g1,g2,cutoff, identity,countDict):
    samplelisttotal=g1 + g2
    import glob,re
    if str(identity) in 'All':
       directory="clust"+str(identity) 
       countsvec={}#initiate counts vector
       countsvec[directory]={} #initiate counts vector
       
       clustervec=getTIGRFAMstruct(['TIGRFAM/corn.fa.hmmout','TIGRFAM/prairie.fa.hmmout'])
       m=0
       for domain in clustervec.keys():
           samplelist=[]
           n=0
           for value in clustervec[domain]:
               contigID=re.sub('_[0-9]$','',value) #remove trailing number
               crow=countDict[contigID+'_'+domain]
               for sample in crow:
                   if sample not in samplelist and countDict[contigID+'_'+domain][sample]>0:
                      samplelist.append(sample)
                      n+=1
           if n>=cutoff: #if the cluster is large enough, get the counts from results file.
              key=1
              countsvec[directory][domain]={}
              countsvec[directory][domain][key]={}
              m+=1
              for value in clustervec[domain]:
                  contigID=re.sub('_[0-9]$','',value) #remove trailing number
                  crow=countDict[contigID+'_'+domain]
                  for sample in crow:
                      if sample in countsvec[directory][domain][key].keys():
                         countsvec[directory][domain][key][sample]+=countDict[contigID+'_'+domain][sample]
                      else:
                           countsvec[directory][domain][key][sample]=countDict[contigID+'_'+domain][sample]
               #check mean count
              meancount=sum(countsvec[directory][domain][key].values())/float(len(samplelist))
              if meancount<10: #remove if mean count is less than 10.
                 del(countsvec[directory][domain])    
       
    else:
         samplelisttotal=g1 + g2
         domainlist=glob.glob('domains/*.fasta') #get a list of all domains
         directory="clust"+str(identity) #loop through each ID cutoff
         print "getting subBins for " + directory
         countsvec={}#initiate counts vector
         countsvec[directory]={} #initiate counts vector
         for domain in domainlist:
             domain=domain.split('/')[1]
             domain=domain.rstrip('.fasta')
             countsvec[directory][domain]={} #initiate counts vector
             clusters=getClusterStruct('domains/'+directory+'/'+domain+'.fasta.uc')
             for key in clusters.keys(): #for each cluster, check if the cluster is large enough for statistical testing
                 #samplelist=[]
                 #n=0
                 #for value in clusters[key]:
                     #sample=value.split('_')[-2]
                     #contigID=re.sub('_[0-9]$','',value) #remove trailing number
                     #crow=countDict[contigID+'_'+domain]
                     #for sample in crow:
                         #if sample not in samplelist:
                            #samplelist.append(sample)
                            #n+=1
                        
                 #if n>=cutoff: #if the cluster is large enough, get the counts from results file.
                 samplelist=[]
                 n=0
                 for value in clusters[key]:
                     contigID=re.sub('_[0-9]$','',value) #remove trailing number
                     crow=countDict[contigID+'_'+domain]
                     for sample in crow:
                         if sample not in samplelist and countDict[contigID+'_'+domain][sample]>0:
                            samplelist.append(sample)
                            n+=1
                 
                 if n>=cutoff:
                    countsvec[directory][domain][key]={}#initiate counts vector
                    totcount=0
                    for value in clusters[key]:
                        contigID=re.sub('_[0-9]$','',value) #remove trailing number
                        crow=countDict[contigID+'_'+domain]
                        for sample in crow:
                            if sample in countsvec[directory][domain][key].keys():
                               countsvec[directory][domain][key][sample]+=countDict[contigID+'_'+domain][sample]
                            else:
                                 countsvec[directory][domain][key][sample]=countDict[contigID+'_'+domain][sample]
                    #check mean count
                    meancount=sum(countsvec[directory][domain][key].values())/float(len(samplelist))
                    if meancount<10: #remove if mean count is less than 10.
                       del(countsvec[directory][domain][key])
               
    
    #countsvec[directory][domain][key][sample]=samplecount
    
    print "writing counts to file for " + directory
    g=open('results_filtered2/abundance_matrix_all_'+directory+'.txt','w')
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

def readCountFile(fileName):
    '''read the tab separated counts file and save into a dictionary'''
    print "reading count file"
    f=open(fileName)
    countDict={}
    for line in f:
        line=line.rstrip()
        line=line.split('\t')
        contigid=line[0]
        contigcount=line[1]
        if contigid not in countDict:
           countDict[contigid]=int(contigcount)
    return countDict

    
def main(fileName,identity):
    #function for calculating the many sub-bins with enough individuals
    #The subbin should contain at least 75% of all individuals. 
    import os,cPickle
    if not os.path.isfile(fileName):
       countDict=getCountStruct('results/')
       f=open(fileName,'wb')
       cPickle.dump(countDict,f)
       f.close()
    else:
         f=open(fileName,'rb')
         countDict=cPickle.load(f)
         f.close()
    
    
    g1=['DLM002','DLM003','DLM004','DLM005','DLM007','DLM008','DLM011','DLM014','DOM005','DOM012','DOM017','DOM020','DOM023','DOM024','DOM025']
    g2=['NLM001','NLM005','NLM021','NLM022','NLM026','NLM029','NLM031','NLM032','NOM001','NOM008','NOM013','NOM017','NOM026','NOM027','NOM028']
    cutoff=round((len(g1)+len(g2))*0.75)
    #read results files
    #countDict=readCountFile(fileName)
    getSubBins(g1,g2,cutoff,identity,countDict)
    
    
    
    

    
if __name__=="__main__":
   import sys
   resultsFileName=sys.argv[1]
   identity=sys.argv[2]
   main(resultsFileName,identity)