#!/bin/bash
#program for clustering, usage 
path=$1
cluster_coeff=$2
filelist=($( ls $path/*.fasta ))
for file in "${filelist[@]}"
    do
    if [ ! -f clust$cluster_coeff/$file.uc ];
    then 
         usearch -cluster_fast $file -id $cluster_coeff -uc clust$cluster_coeff/$file.uc
    fi
done
      