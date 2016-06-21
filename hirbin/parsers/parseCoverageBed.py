#!/usr/bin/env python
def parseCoverageBed(filename,outputfilename):
  f=open(filename)
  g=open(outputfilename,'w')
  for line in f:
    line=line.rstrip()
    line=line.split('\t')
    if len(line)>7: #extra check
      ref=line[0]
      start=line[3]
      stop=line[4]
      strand=line[6]
      idx=line[8].lstrip('ID=')
      counts=line[-1]
      newline=ref+'_'+idx+':'+start+':'+stop+':'+strand+'\t'+counts+'\n'
      g.write(newline)
  f.close()
  g.close()

if __name__=='__main__':
  import sys
  parseCoverageBed(sys.argv[1],sys.argv[2])