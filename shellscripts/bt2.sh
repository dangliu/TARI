#!/bin/bash

cd /home/tari/tomato/RAD2Core/alignments/20180704   
SAMPLE=$1
INDEX=$(echo $SAMPLE | sed 's/_r[0-9]*//')
LANE=$(echo ${INDEX:3:1})
echo "Working on '$SAMPLE'."
fq1=../../cleaned/20180704/$INDEX/$SAMPLE.1.fq.gz  
fq2=../../cleaned/20180704/$INDEX/$SAMPLE.2.fq.gz  
bam_file=$INDEX/$SAMPLE.bam
bt2_db=../../genome/ref
# PU needs to be modified alone with different sequencing runs
# -p can be modified by using variant number of threads
bowtie2 -x $bt2_db -p 8 --very-sensitive --dovetail --no-mixed --no-discordant -X 800 --rg PU:CCK2DANXX.$LANE --rg PG:bowtie2 --rg PL:ILLUMINA --rg SM:$SAMPLE --rg-id $SAMPLE -1 $fq1 -2 $fq2 | samtools view -b  -o $bam_file
