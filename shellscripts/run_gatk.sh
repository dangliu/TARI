#!/bin/bash

# sam to sorted.bam
#for i in *bam
#do
#	x=$(echo $i | sed 's/\.bam//')
#	echo $x.sort
#	gatk SortSam \
#	-I $i \
#	-O $x.sorted.bam \
#	-SO coordinate 
#done


## markduplicate
#for i in *.sorted.bam
#do
#	x=$(echo $i | sed 's/\.sorted\.bam//')
#	echo $x.dedup
#	gatk MarkDuplicates \
#	-I $i \
#	-O $x.dedup.bam \
#	-M $x.metrics.txt \
#	--TAGGING_POLICY All	
#done

#for i in *.dedup.bam
#do 
#	bamtools stats -in $i > $i.stats
#done

## Make a query name list of reads to be included
#for i in *dedup.bam
#do
#	x=$(echo $i | sed 's/\.dedup\.bam//')
#	echo $x.includeReadList
#	Rscript ~/bin/Rscript/FilterBAM.R $i $x.IncludeQueryName.txt
#done

## Filter bam according to IncludeQueryName
#for i in *dedup.bam
#do
#	x=$(echo $i | sed 's/\.dedup\.bam//')
#	echo $x.includeReadList	
#	gatk FilterSamReads \
#	--FILTER includeReadList \
#	-I $i \
#	-O $x.filtered.bam \
#	-RLF $x.IncludeQueryName.txt \
#	-SO coordinate 
#done

## Index bam
#for i in *.filtered.bam
#do 
#	samtools index $i
#done

## Call haplotype
#for i in *.filtered.bam
#do
#	x=$(echo $i | sed 's/\.filtered\.bam//')
#	echo $x.HC
#	gatk HaplotypeCaller \
#	-R ~/tomato/RAD2Core/genome/ref.normalized.fa \
#	-I $i \
#	-O $x.g.vcf.gz \
#	-ERC GVCF 
#done


### For making dbsnp and finall snp set

## Load in GenomicDB and join-genotype
#echo Import GenomicDBs and Genotype cohort GVCFs...
#for i in $(seq 1 12)
#do
#	if [ $i -lt 10 ]
#	then
#		echo processing SL3.0ch0$i
#		gatk GenomicsDBImport --sample-name-map ../info/QF01.map  --genomicsdb-workspace-path chr0$i.database -L SL3.0ch0$i
#		gatk GenotypeGVCFs -R ~/tomato/RAD2Core/genome/ref.normalized.fa -V gendb://chr0$i.database -O chr0$i.vcf.gz
#	else
#		echo processing SL3.0ch$i
#		gatk GenomicsDBImport --sample-name-map ../info/QF01.map  --genomicsdb-workspace-path chr$i.database -L SL3.0ch$i
#		gatk GenotypeGVCFs -R ~/tomato/RAD2Core/genome/ref.normalized.fa -V gendb://chr$i.database -O chr$i.vcf.gz
#	fi
#done


## Combine VCF
# echo Gather VCFs...
# gatk GatherVcfs \
# -I chr01.vcf.gz \
# -I chr02.vcf.gz \
# -I chr03.vcf.gz \
# -I chr04.vcf.gz \
# -I chr05.vcf.gz \
# -I chr06.vcf.gz \
# -I chr07.vcf.gz \
# -I chr08.vcf.gz \
# -I chr09.vcf.gz \
# -I chr10.vcf.gz \
# -I chr11.vcf.gz \
# -I chr12.vcf.gz \
# -O combined.vcf.gz

## Index vcf
# tabix -p vcf combined.vcf.gz

## Filter variants
# echo Vilter variants...
# gatk VariantFiltration \
# -R ~/tomato/RAD2Core/genome/ref.normalized.fa \
# -V combined.vcf.gz \
# -filter "QUAL < 1000.0 || QD < 20.0 || FS > 60.0 || SOR > 3.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || AN < 768" \
# --filter-name "my_snp_filter" \
# -O combined.filtered.vcf.gz

## Select SNP variants
# echo Select variants...
# gatk SelectVariants \
# -R ~/tomato/RAD2Core/genome/ref.normalized.fa \
# -V combined.filtered.vcf.gz \
# -select-type SNP \
# --restrict-alleles-to BIALLELIC \
# --exclude-filtered \
# -O combined.filtered.snp.vcf.gz

## Select INDEL variants
# gatk SelectVariants -R ~/tomato/RAD2Core/genome/ref.normalized.fa -V combined.filtered.vcf.gz -select-type INDEL --exclude-filtered -O combined.filtered.indel.vcf.gz




## Base recalibaration
# for i in *filtered.bam
# do
#         x=$(echo $i | sed 's/\.filtered\.bam//')
#         echo $x.recal
#     	gatk BaseRecalibrator \
#         -R ~/tomato/RAD2Core/genome/ref.normalized.fa \
#         -I $i \
#         --known-sites ~/tomato/RAD2Core/gatk/info/db.snp.vcf.gz \
#         --known-sites ~/tomato/RAD2Core/gatk/info/db.indel.vcf.gz \
#         -O $x.recal.table
# done

# for i in *filtered.bam
# do
#         x=$(echo $i | sed 's/\.filtered\.bam//')
#         echo $x.applyBQSR
# 		gatk ApplyBQSR \
#         -R ~/tomato/RAD2Core/genome/ref.normalized.fa \
#         -I $i \
#         -bqsr $x.recal.table \
#         -O $x.recal.bam
# done


## 2nd haplotype-calling
for i in *.recal.bam
do
	x=$(echo $i | sed 's/\.recal\.bam//')
	echo $x.2nd.HC
	gatk HaplotypeCaller \
	-R ~/tomato/RAD2Core/genome/ref.normalized.fa \
	-I $i \
	-O $x.2nd.g.vcf.gz \
	-ERC GVCF 
done
