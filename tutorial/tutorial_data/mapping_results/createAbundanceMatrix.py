#!/usr/bin/env python

import sys
import glob
def parseargs():
  ''' Function for parsing the argument vector'''
  args=sys.argv[1:]
  if len(args)==1:
    #argument is a pattern or only one file
    filelist=glob.glob(args[0])
  else:
    #argument is several filenames separated by spaces
    filelist=args
  return(filelist)

def main(filelist):
  ''' Function for parsing the output from the mapping in TIGRFAM format and creating an abundance matrix'''
  
  #find unique samples
  samplelist=[]
  domains={}
  #create sample names
  for filename in filelist:
      sample=filename.split('/')[-1].split('.')[0]
      if sample not in samplelist:
         samplelist.append(sample)
  for filename in filelist: #loop through files and pick out TIGRFAMs + counts
      sample=filename.split('/')[-1].split('.')[0]
      f=open(filename)
      for line in f:
          name=line.split('\t')[0]
          count=int(line.split('\t')[1])
          tigrfam=name.split('_')[2].split(':')[0]
          if tigrfam in str(range(0,10)): #check if TIGRFAM name ends with _1
            tigrfam=name.split('_')[3].split(':')[0]
          if tigrfam not in domains:
             domains[tigrfam]={}
             for i in samplelist:
                 domains[tigrfam][i]=0
          domains[tigrfam][sample] += count
      f.close()
  #save results in the file abundance_matrix.txt
  g=open('abundance_matrix.txt','w')
  g.write('domain')
  for sample in samplelist:
      g.write('\t'+sample)
  g.write('\n')
  for domain in domains.keys():
      g.write(domain)
      for sample in samplelist:
          g.write('\t'+str(domains[domain][sample]))
      g.write('\n')
  g.close()
    
                                  

if __name__=='__main__':
  filelist=parseargs()
  main(filelist)