
def main(referencedir):
    import glob,gzip
    filelist=glob.glob(referencedir+'*.fasta')
    from Bio import SeqIO
    domains={}
    for filename in filelist:
        f1=open(filename)
        sequences = SeqIO.parse(f1, "fasta")
        for s in sequences:
            domain=s.description.split()[1]
            if domain not in domains.keys():
               domains[domain]=[]
               domains[domain].append(s)
            else:
                 domains[domain].append(s)
        f1.close()
    
    for domain in domains:
        g=open('domains/'+domain+'.fasta','w')
        SeqIO.write(domains[domain], g, "fasta")
        g.close()
    
if __name__=='__main__':
   import sys
   referencedir=sys.argv[1]
   main(referencedir)