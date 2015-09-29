#!/bin/bash
samplelist=(
NOM017
NOM001
NOM027
NOM028

)
for i in "${samplelist[@]}"
do
  python convertCoord2.py $i
done

