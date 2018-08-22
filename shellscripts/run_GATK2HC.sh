#!/bin/bash

# Sort sam and convert to sorted.bam
for i in *bam
do
	x=$(echo $i | sed 's/\.bam//')
	echo $x.sort
	gatk SortSam \
	-I $i \
	-O $x.sorted.bam \
	-SO coordinate 
done

## Markduplicate
for i in *.sorted.bam
do
	x=$(echo $i | sed 's/\.sorted\.bam//')
	echo $x.dedup
	gatk MarkDuplicates \
	-I $i \
	-O $x.dedup.bam \
	-M $x.metrics.txt \
	--TAGGING_POLICY All	
done

## Make a stats after de-duplication
for i in *.dedup.bam
do 
	bamtools stats -in $i > $i.stats
done

## Make a query name list of reads to be included
for i in *dedup.bam
do
	x=$(echo $i | sed 's/\.dedup\.bam//')
	echo $x.includeReadList
	Rscript ~/bin/Rscript/FilterBAM.v2.R $i $x.IncludeQueryName.txt
done

## Filter bam according to IncludeQueryName
for i in *dedup.bam
do
	x=$(echo $i | sed 's/\.dedup\.bam//')
	echo $x.includeReadList	
	gatk FilterSamReads \
	--FILTER includeReadList \
	-I $i \
	-O $x.filtered.bam \
	-RLF $x.IncludeQueryName.txt \
	-SO coordinate 
done

## Make a stats after bam-filtering
for i in *.filtered.bam
do 
	bamtools stats -in $i > $i.stats
done

## Base recalibaration
for i in *filtered.bam
do
        x=$(echo $i | sed 's/\.filtered\.bam//')
        echo $x.recal
    	gatk BaseRecalibrator \
        -R ~/tomato/RAD2Core/genome/ref.normalized.fa \
        -I $i \
        --known-sites ~/tomato/RAD2Core/gatk/info/db.snp.vcf.gz \
        --known-sites ~/tomato/RAD2Core/gatk/info/db.indel.vcf.gz \
        -O $x.recal.table
done

## Apply BQSR on filtered BAM
for i in *filtered.bam
do
        x=$(echo $i | sed 's/\.filtered\.bam//')
        echo $x.applyBQSR
		gatk ApplyBQSR \
        -R ~/tomato/RAD2Core/genome/ref.normalized.fa \
        -I $i \
        -bqsr $x.recal.table \
        -O $x.recal.bam
done

## Call haplotype
for i in *.recal.bam
do
	x=$(echo $i | sed 's/\.recal\.bam//')
	echo $x.HC
	gatk HaplotypeCaller \
	-R ~/tomato/RAD2Core/genome/ref.normalized.fa \
	-I $i \
	-O $x.g.vcf.gz \
	-ERC GVCF 
done
