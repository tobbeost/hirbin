#!/usr/bin/env python2.7
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
                if target not in clusters[clusterID]:
                   clusters[clusterID].append(target)
       
    f.close()
    return clusters

def getTIGRFAMstruct(domain):
    filename='domains/'+domain+'.fasta'
    f=open(filename)
    clusters=[]
    for line in f:
        if line.startswith('>'):
           line=line.lstrip('>').split()
           line=line[0]
           if line not in clusters:
              clusters.append(line)
    f.close()
    return(clusters)
    
def getContigs(clustdict,domain):
    countsvec={}#initiate counts vector
    countsvec['All']=getTIGRFAMstruct(domain)
    for identity in clustdict:
        directory="clust"+str(identity) #loop through each ID cutoff  
        countsvec[directory]=[] #initiate counts vector
        clusters=getClusterStruct('domains/'+directory+'/'+domain+'.fasta.uc')
        countsvec[identity]=clusters[int(clustdict[identity])]
    return(countsvec)
    print "done"

def readFile(fileName):
    f=open(fileName)
    clusters={}
    for line in f:
        line=line.rstrip()
        line=line.split()
        if len(line)>1:
           clusters[line[0]]=line[1]
    f.close()
    return(clusters)
    
def checkConsistancy(contigs,clusters,domain):
    import collections
    clustlist=[0.95,0.9,0.85,0.8,0.75,0.7,0.65,0.6,0.55,0.5,0.45,0.4,'All']
    clustlist2=[str(x) for x in clustlist]
    ids=[]
    clustlist2.reverse()
    for identity in clustlist2:
        if identity in contigs:
           ids.append(identity)
    #g=open('cluster_overlap_'+domain+'_5.txt','w')
    print '\t'.join(['identity','cluster','size','overlap from previous','new from previous','disappeared from previous'])+'\n',
    print 'All\tNA\t'+str(len(contigs['All']))+'\t-\t-\t-\n',
    t1=ids[0:(len(ids)-1)]
    t2=ids[1:len(ids)]
    for a,b in zip(t1,t2):
        alist=collections.Counter(contigs[a])
        blist=collections.Counter(contigs[b])
        overlap = list((alist & blist))
        a_remainder = list((alist - blist))
        b_remainder = list((blist - alist))
        print b+'\t'+clusters[b]+'\t'+str(len(contigs[b]))+'\t'+str(len(overlap))+'\t'+str(len(b_remainder))+'\t'+str(len(a_remainder))+'\n',
        
    #g.close()
           
def main(fileName,domain):
    clusters=readFile(fileName)
    contigs=getContigs(clusters,domain)
    checkConsistancy(contigs,clusters,domain)
    
    
  
    
    

    
if __name__=="__main__":
   import sys
   fileName=sys.argv[1]
   domain=sys.argv[2]
   main(fileName,domain)