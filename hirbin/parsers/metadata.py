# coding: utf-8
"""

"""
import os

class Hirbin_run(object):
  def __init__(self,name):
    self.name=name
      
  def getName(self):
    return self.name

  def createOutputDirectory(self, output_directory):
    ''' Function for creating a new output directory. If the name is not specified, decide a new name '''
    if output_directory==None: #no output directory specified, create a new output directory
      new_output_dir='hirbin_output'
      not_created_dir=False
      suffix=1
      #create a new name for the output directory
      while not not_created_dir:
        if os.path.isdir(new_output_dir):
          suffix=suffix+1
          new_output_dir='hirbin_output'+str(suffix)
        else:
          not_created_dir=True
      try:
        os.mkdir(new_output_dir)
        output_directory=new_output_dir
        print "Creating a new output directory at "+os.getcwd()+'/'+output_directory
      except OSError as e:
        if not force:
          raise
    else:
      if not os.path.isdir(output_directory):
        try:
          os.mkdir(output_directory)
        except OSError as e:
          print "Output directory already exists, you can use an already existing output directory by including the flag -f"
          raise
    #print output_directory
    self.output_directory=output_directory
    return(output_directory)
  
  def readMetadata(self,metadataFile):
    '''
    Reads in the sample mapping/metadata file and creates a dictionary 
    with sample names and file path to assembly for that sample.
    '''
    with open(metadataFile) as f:
      line=f.readline() #read header
      line=line.lower() #ignore case
      header=line.split()
      if 'name' in header:
        nameindex=header.index('name')
        self.samples=[]
      else:
        self.samples=None
      if 'group' in header:
        groupsindex=header.index('group')
        self.groups={}
      else:
        self.groups=None
      if 'reference' in header:
        referenceindex=header.index('reference')
        self.reference={}
      else:
        self.reference=None
      if 'annotation' in header:
        annotationindex=header.index('annotation')
        self.annotation={}
      else:
        self.annotation=None
      if 'counts' in header:
        countsindex=header.index('counts')
        self.counts={}
      else:
        self.counts=None
      if 'reads1' in header:
        reads1index=header.index('reads1')
        self.reads1={}
      else:
        self.reads1=None
      if 'reads2' in header:
        reads2index=header.index('reads2')
        self.reads2={}
      else:
        self.reads2=None
        
      for line in f:
        line=line.rstrip()
        line=line.split()
        samplename=line[nameindex]
        if not (self.samples is None):
          if samplename not in self.samples:
            self.samples.append(samplename)
        if not (self.groups is None):
          gr=line[groupsindex]
          if samplename not in self.groups:
            self.groups[samplename]=gr
        if not (self.reference is None):
          ref=line[referenceindex]
          if samplename not in self.reference:
            self.reference[samplename]=ref
        if not (self.annotation is None):
          an=line[annotationindex]
          if samplename not in self.annotation:
            self.annotation[samplename]=an
        if not (self.counts is None):
          co=line[countsindex]
          if samplename not in self.counts:
            self.counts[samplename]=co
        if not (self.reads1 is None):
          re1=line[reads1index]
          if samplename not in self.reads1:
            self.reads1[samplename]=re1
        if not (self.reads2 is None):
          re2=line[reads2index]
          if samplename not in self.reads2:
            self.reads2[samplename]=re2

    def getSamples(self):
      return self.groups.keys()
    def getGroups(self):
      return self.groups
    def getAnnotation(self):
      return self.annotation
    def getReference(self):
      return self.reference
    def getCounts(self):
      return self.counts
    def getReads1(self):
      return self.reads1
    def getReads2(self):
      return self.reads2