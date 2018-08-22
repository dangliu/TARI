#!/bin/bash
cd /home/tari/tomato/RAD2Core/stacks
mkdir BAM.wfilter/$1/
ln -s ~/home/tari/tomato/RAD2Core/gatk/BAM/$1/*.filtered.bam ./BAM.wfilter/$1/
rename 's/\.filtered//' BAM.wfilter/$1/*
mkdir gstaks.wfilter/$1/
cat popmap.txt | sed 's/QF0101/'$1'/g' >$1.popmap.txt
gstacks -I BAM.wfilter/$1/ -O gstaks.wfilter/$1/ -M $1.popmap.txt -t 16 --rm-pcr-duplicates >gstaks.wfilter/$1/$1.gstacks.wfilter.oe
cd gstaks.wfilter/$1/
populations -P ./ -M ../../$1.popmap.txt -r 0.80 --min_maf 0.05 --max_obs_het 0.7 --vcf --genepop -t 16 >$1.populations.oe
rm ../../$1.popmap.txt
