#!/bin/bash

for i in *g.vcf.gz
do
	x=$(echo $i | sed 's/\.g\.vcf\.gz//')
	echo $x	
	gatk ValidateVariants \
	-R ~/tomato/RAD2Core/genome/ref.normalized.fa \
	-V $i 
done
