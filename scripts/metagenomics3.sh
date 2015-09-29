#!/bin/bash
seqlib=$1
evaluecutoff=$2
filelist=($(ls $seqlib *.fastq))
echo "${#filelist[@]}"
echo $numthreads
if [ ! -d "TIGRFAM" ]; then
  mkdir TIGRFAM
fi
for i in "${filelist[@]}"
do
    echo $seqlib/$i
    if [ ! -f TIGRFAM/$i.hmmout ];
    then
        zcat $seqlib/$i | transeq --frame 6 --filter > $seqlib/protseq/$i.pep
        #run TIGRFAM
        hmmsearch --domtblout TIGRFAM/$i.hmmout -E $evaluecutoff /storage/shared/db/TIGRFAM/13.0/TIGRFAM.HMM $seqlib/protseq/$i.pep > tmpout.out
    fi
done
echo "done"