
class Hirbin_run(object):
  def __init__(self,name):
    self.name=name
  def getName(self):
    return self.name
  
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
        
      samples={}
      for line in f:
        line=line.rstrip()
        line=line.split()
        samplename=line[nameindex]
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

    def getSamples(self):
      return self.groups.keys()
    def getGroups(self):
      return self.groups
    def getAnnotation(self):
      return self.annotation
    def getReference(self):
      return self.reference