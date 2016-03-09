#!/usr/bin/env python
from extractSequences2 import *
import random

def randomdomains(filename):
  f=open(filename)
  domainlist=[]
  for line in f:
    if line.startswith('TIGR'):
      line=line.split('\t')[0]
      domainlist.append(line)
  n=len(domainlist)
  newn=round(n*0.01)
  newdomainlist=random.sample(domainlist,int(newn))
  g=open('selected_TIGRFAMS.txt','w')
  for domain in newdomainlist:
    g.write(domain+'\n')
  f.close()
  g.close()
    
def getSelectedDomains(filename):
  f=open(filename)
  domainlist=[]
  for line in f:
    line=line.rstrip()
    domainlist.append(line)
  return(domainlist)

def reduceannot(sample,domainlist):
  f=open('TIGRFAM/'+sample+'.contigs.m100.gz.hmmout')
  g=open('TIGRFAM/reduced/'+sample+'.hmmout','w')
  f.readline()
  for line in f:
    parts=line.split()
    if len(parts)>4 and parts[3] in domainlist:
      g.write(line)
  f.close()
  g.close()

def reducefastafile(oldfilename,newfilename,contiglist):
  f=open(oldfilename)
  g=open(newfilename,'w')
  printseq=False
  for line in f:
    if line.startswith('>'):
      printseq=False
      idx=line.lstrip('>').split()[0]
      if idx in contiglist:
        printseq=True
        g.write(line)
    elif printseq==True:
      g.write(line)
  f.close()
  g.close()
  

def reducecontigs(sample):
  (annot,start,stop,strand,evalue)=load_annotation_pfam('TIGRFAM/reduced/'+sample+'.hmmout')
  protlist=[]
  nucllist=[]
  ntot=0
  for domain in annot:
    for contigprot in annot[domain]:
      ntot+=1
      if contigprot not in protlist:
        protlist.append(contigprot)
      nuclseq='_'.join(contigprot.split('_')[0:2])
      if nuclseq not in nucllist:
        nucllist.append(nuclseq)
  #print protlist[0:10],len(protlist)
  #print nucllist[0:10],len(nucllist)
  #print ntot
  reducefastafile('contigs/'+sample+'.contigs.m100','contigs/reduced/'+sample+'.fasta',nucllist)
  reducefastafile('protseq/'+sample+'.contigs.m100.gz.pep','protseq/reduced/'+sample+'.fasta',protlist)
  return(nucllist)

def reducemappingresults(sample,nucllist,domainlist):
  f=open('mapping_results/'+sample+'.fastq.tab')
  g=open('mapping_results/reduced/'+sample+'.fastq.tab','w')
  for line in f:
    idd=line.split(':')[0]
    parts=idd.split('_')
    contig='_'.join(parts[0:2])
    domain=parts[-1]
    if contig in nucllist and domain in domainlist:
      g.write(line)
  f.close()
  g.close()

def main():
  #randomdomains('mapping_results/abundance_matrix.txt')
  domainlist=getSelectedDomains('selected_TIGRFAMS.txt')
  #print len(domainlist)
  samplelist=['DLM005','DOM025','DOM012','DOM005','DOM017','NLM021','NLM031','NOM017','NOM008','NOM026']
  #for sample in samplelist:
  #  reduceannot(sample,domainlist)
  sample='DLM005'
  for sample in samplelist:
    nucllist=reducecontigs(sample)
    reducemappingresults(sample,nucllist,domainlist)
  
    
  #(annot,start,stop,strand,evalue)=load_annotation_pfam('TIGRFAM/reduced/'+sample+'.hmmout')
  #print len(annot.keys())
  #print annot.keys()[0]
  
  #print annot[annot.keys()[0]]

if __name__=='__main__':
  main()