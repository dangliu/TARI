#!/bin/bash


## Load in GenomicDB and join-genotype
# nohup gatk GenomicsDBImport --batch-size 50 --sample-name-map ../info/all.v2.map --genomicsdb-workspace-path chr01.database -L SL3.0ch01 >chr01.GDB.oe &
# nohup gatk GenomicsDBImport --batch-size 50 --sample-name-map ../info/all.v2.map --genomicsdb-workspace-path chr02.database -L SL3.0ch02 >chr02.GDB.oe &
# nohup gatk GenomicsDBImport --batch-size 50 --sample-name-map ../info/all.v2.map --genomicsdb-workspace-path chr03.database -L SL3.0ch03 >chr03.GDB.oe &
# nohup gatk GenomicsDBImport --batch-size 50 --sample-name-map ../info/all.v2.map --genomicsdb-workspace-path chr04.database -L SL3.0ch04 >chr04.GDB.oe &
# nohup gatk GenomicsDBImport --batch-size 50 --sample-name-map ../info/all.v2.map --genomicsdb-workspace-path chr05.database -L SL3.0ch05 >chr05.GDB.oe &
# nohup gatk GenomicsDBImport --batch-size 50 --sample-name-map ../info/all.v2.map --genomicsdb-workspace-path chr06.database -L SL3.0ch06 >chr06.GDB.oe &
# nohup gatk GenomicsDBImport --batch-size 50 --sample-name-map ../info/all.v2.map --genomicsdb-workspace-path chr07.database -L SL3.0ch07 >chr07.GDB.oe &
# nohup gatk GenomicsDBImport --batch-size 50 --sample-name-map ../info/all.v2.map --genomicsdb-workspace-path chr08.database -L SL3.0ch08 >chr08.GDB.oe &
# nohup gatk GenomicsDBImport --batch-size 50 --sample-name-map ../info/all.v2.map --genomicsdb-workspace-path chr09.database -L SL3.0ch09 >chr09.GDB.oe &
# nohup gatk GenomicsDBImport --batch-size 50 --sample-name-map ../info/all.v2.map --genomicsdb-workspace-path chr10.database -L SL3.0ch10 >chr10.GDB.oe &
# nohup gatk GenomicsDBImport --batch-size 50 --sample-name-map ../info/all.v2.map --genomicsdb-workspace-path chr11.database -L SL3.0ch11 >chr11.GDB.oe &
# nohup gatk GenomicsDBImport --batch-size 50 --sample-name-map ../info/all.v2.map --genomicsdb-workspace-path chr12.database -L SL3.0ch12 >chr12.GDB.oe &

# nohup gatk GenotypeGVCFs -R ~/tomato/RAD2Core/genome/ref.normalized.fa -V gendb://chr01.database -O chr01.vcf.gz > chr01.GT.oe &
# nohup gatk GenotypeGVCFs -R ~/tomato/RAD2Core/genome/ref.normalized.fa -V gendb://chr02.database -O chr02.vcf.gz > chr02.GT.oe &
# nohup gatk GenotypeGVCFs -R ~/tomato/RAD2Core/genome/ref.normalized.fa -V gendb://chr03.database -O chr03.vcf.gz > chr03.GT.oe &
# nohup gatk GenotypeGVCFs -R ~/tomato/RAD2Core/genome/ref.normalized.fa -V gendb://chr04.database -O chr04.vcf.gz > chr04.GT.oe &
# nohup gatk GenotypeGVCFs -R ~/tomato/RAD2Core/genome/ref.normalized.fa -V gendb://chr05.database -O chr05.vcf.gz > chr05.GT.oe &
# nohup gatk GenotypeGVCFs -R ~/tomato/RAD2Core/genome/ref.normalized.fa -V gendb://chr06.database -O chr06.vcf.gz > chr06.GT.oe &
# nohup gatk GenotypeGVCFs -R ~/tomato/RAD2Core/genome/ref.normalized.fa -V gendb://chr07.database -O chr07.vcf.gz > chr07.GT.oe &
# nohup gatk GenotypeGVCFs -R ~/tomato/RAD2Core/genome/ref.normalized.fa -V gendb://chr08.database -O chr08.vcf.gz > chr08.GT.oe &
# nohup gatk GenotypeGVCFs -R ~/tomato/RAD2Core/genome/ref.normalized.fa -V gendb://chr09.database -O chr09.vcf.gz > chr09.GT.oe &
# nohup gatk GenotypeGVCFs -R ~/tomato/RAD2Core/genome/ref.normalized.fa -V gendb://chr10.database -O chr10.vcf.gz > chr10.GT.oe &
# nohup gatk GenotypeGVCFs -R ~/tomato/RAD2Core/genome/ref.normalized.fa -V gendb://chr11.database -O chr11.vcf.gz > chr11.GT.oe &
# nohup gatk GenotypeGVCFs -R ~/tomato/RAD2Core/genome/ref.normalized.fa -V gendb://chr12.database -O chr12.vcf.gz > chr12.GT.oe &



## Combine VCF
echo Gather VCFs...
gatk GatherVcfs \
-I chr01.vcf.gz \
-I chr02.vcf.gz \
-I chr03.vcf.gz \
-I chr04.vcf.gz \
-I chr05.vcf.gz \
-I chr06.vcf.gz \
-I chr07.vcf.gz \
-I chr08.vcf.gz \
-I chr09.vcf.gz \
-I chr10.vcf.gz \
-I chr11.vcf.gz \
-I chr12.vcf.gz \
-O combined.vcf.gz

## Index vcf
tabix -p vcf combined.vcf.gz

## Filter variants
echo Vilter variants...
gatk VariantFiltration \
-R ~/tomato/RAD2Core/genome/ref.normalized.fa \
-V combined.vcf.gz \
-filter "QD < 2.0 || FS > 25.0 || SOR > 7.0 || MQ < 20.0 || MQRankSum < -5.0 || ReadPosRankSum < -5.0" \
--filter-name "my_snp_filter" \
-O combined.filtered.vcf.gz

## Select SNP variants
echo Select variants...
gatk SelectVariants \
-R ~/tomato/RAD2Core/genome/ref.normalized.fa \
-V combined.filtered.vcf.gz \
-select-type SNP \
--restrict-alleles-to BIALLELIC \
--exclude-filtered \
-O combined.filtered.snp.vcf.gz

## Select INDEL variants
gatk SelectVariants -R ~/tomato/RAD2Core/genome/ref.normalized.fa -V combined.filtered.vcf.gz -select-type INDEL --exclude-filtered -O combined.filtered.indel.vcf.gz