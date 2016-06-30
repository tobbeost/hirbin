#!/usr/bin/env python
# coding: utf-8
def getCountStruct(metadata):
    '''
    function for reading all results from mapping and save in a dictonary used for clustering
    '''
    print "Reading and saving mapping results..."
    filelist=metadata.counts
    countDict={}

    for sample in filelist:
        f=open(filelist[sample])
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

def writeAbundanceMatrix(metadata,domains):
  '''Function for writing the abundance matrix (for bins) to a file in the hirbin output directory'''
  #save results in the file abundance_matrix.txt
  filelist=metadata.counts
  samplelist=metadata.samples
  output_directory=metadata.output_directory
  g=open(output_directory+'/abundance_matrix_bins.txt','w')
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

def filterAbundanceMatrix(metadata,domains,p,minMeanCount):
  '''Function for filtering the abundance matrix (for bins) based on minimum coverage and representation.'''
  samplelist=metadata.samples
  cutoff=round(len(samplelist)*p)
  filtereddomains={}
  for tigrfam in domains:
    meancount=sum(domains[tigrfam].values())/float(len(samplelist))
    if meancount>=minMeanCount: #remove if mean count is less than minMeanCount.
      n=0
      for sample in domains[tigrfam]:
        if domains[tigrfam][sample]>0:
          n+=1
      if n>=cutoff:
        #passed criteria
        filtereddomains[tigrfam]={}
        for sample in samplelist:  
          filtereddomains[tigrfam][sample]=domains[tigrfam][sample]
  return filtereddomains
      

def createAbundanceMatrix(metadata,p,minMeanCount):
  ''' Function for parsing the output from the mapping in TIGRFAM format and creating an abundance matrix'''
  
  filelist=metadata.counts
  samplelist=metadata.samples
  domains={}
  #create sample names
  for sample in filelist: #loop through files and pick out TIGRFAMs + counts
      f=open(filelist[sample])
      for line in f:
          name=line.split('\t')[0]
          count=int(line.split('\t')[1])
          tigrfam=name.split(':')[0].split('_')[-1]
          if tigrfam in str(range(0,10)): #check if TIGRFAM name ends with _1
            tigrfam=name.split(':')[0].split('_')[-2]
          if tigrfam not in domains:
             domains[tigrfam]={}
             for i in samplelist:
                 domains[tigrfam][i]=0
          domains[tigrfam][sample] += count
      f.close()
  filtereddomains=filterAbundanceMatrix(metadata,domains,p,minMeanCount)
  writeAbundanceMatrix(metadata,filtereddomains)
  return filtereddomains

