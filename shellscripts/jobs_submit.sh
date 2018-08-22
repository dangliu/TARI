#!/bin/bash

cd /home/tari/tomato/RAD2Core/alignments/20180704

# e.g. QF0101
INDEX=$1


# used to submit jobs for bt2.sh

for spl in $(cut -f1 ./$INDEX/sample.txt)
do 
./bt2.sh $spl >>./$INDEX/log 2>>./$INDEX/log2
done
