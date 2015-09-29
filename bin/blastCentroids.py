def getCentroidContigID(domain,clusterid,identity):
    identity=identity.rstrip('0')
    filename='domains/clust'+identity+'/'+domain+'.fasta.uc'
    f=open(filename)
    found=False
    for line in f:
        if line.startswith('S'):
           line=line.split('\t')
           if line[1] in clusterid:
              contigID=line[8]
              found=True
              
    if found==False:
       print 'cluster ' + clusterid + ' not found in file ' + filename
       return('')
    else:
         return(contigID)

def createFastaFile(contigID,domain):
    from Bio import SeqIO
    filename='domains/'+domain+'.fasta'
    f=open(filename)
    candidates=[]
    g=open('tmp.fasta','w')
    for record in SeqIO.parse(f,'fasta'):
        if record.id in contigID:
           candidates.append(record)
    f.close()
    lmax=0
    maxcandidate=''
    for candidate in candidates:
        if len(candidate.seq)>lmax:
           maxcandidate=candidate
    SeqIO.write(maxcandidate, g, "fasta")
    g.close()

def runBlast(outputfilename,nrpath):
    import os
    line='blastp -query tmp.fasta -db '+nrpath+' -evalue 1e-10 -out '+outputfilename+' -num_alignments 30 -num_threads 10 -outfmt "6 qseqid stitle skingdom pident evalue bitscore"'
    print 'running blast and saving to file ' + outputfilename
    os.system(line)
    
def main(domain,clusterid,identity,nrpath):
    contigID=getCentroidContigID(domain,clusterid,identity)
    createFastaFile(contigID,domain)
    outputfilename='blastoutput/'+domain+'_clust'+identity+'_'+clusterid+'.out'
    runBlast(outputfilename,nrpath)


if __name__=='__main__':
   import sys
   args=sys.argv[1:]
   domain=args[0]
   clusterid=args[1]
   identity=args[2]
   nrpath=args[3]
   main(domain,clusterid,identity,nrpath)