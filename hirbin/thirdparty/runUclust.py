#!/usr/bin/env python
# coding: utf-8
import os
def runUclust(path,identityCutoff):
  print "Clustering sequences using sequence identity cutoff " +str(identityCutoff)+'.' 
  filelist=os.listdir(path)
  outpath=path+'../clust'+str(identityCutoff)+'/'
  for f in filelist:
    if not os.path.isfile(outpath+f+'.uc'):
      #if aborted continue at previous stage
      command="usearch -cluster_fast "+ path+f+" -id "+ str(identityCutoff) + " -uc " +outpath+f+'.uc &> /dev/null'
      os.system(command)

if __name__=='__main__':
  import sys
  runUclust(sys.argv[1],sys.argv[2])

      