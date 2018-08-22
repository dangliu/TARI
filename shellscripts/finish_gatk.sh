#!/bin/bash


## Load in GenomicDB and join-genotype
# nohup gatk GenomicsDBImport --batch-size 50 --sample-name-map ../info/all.map --genomicsdb-workspace-path chr01.database -L SL3.0ch01 >chr1.GDB.oe &
# nohup gatk GenomicsDBImport --batch-size 50 --sample-name-map ../info/all.map --genomicsdb-workspace-path chr02.database -L SL3.0ch02 >chr2.GDB.oe &
# nohup gatk GenomicsDBImport --batch-size 50 --sample-name-map ../info/all.map --genomicsdb-workspace-path chr03.database -L SL3.0ch03 >chr3.GDB.oe &
# nohup gatk GenomicsDBImport --batch-size 50 --sample-name-map ../info/all.map --genomicsdb-workspace-path chr04.database -L SL3.0ch04 >chr4.GDB.oe &
# nohup gatk GenomicsDBImport --batch-size 50 --sample-name-map ../info/all.map --genomicsdb-workspace-path chr05.database -L SL3.0ch05 >chr5.GDB.oe &
# nohup gatk GenomicsDBImport --batch-size 50 --sample-name-map ../info/all.map --genomicsdb-workspace-path chr06.database -L SL3.0ch06 >chr6.GDB.oe &
# nohup gatk GenomicsDBImport --batch-size 50 --sample-name-map ../info/all.map --genomicsdb-workspace-path chr07.database -L SL3.0ch07 >chr7.GDB.oe &
# nohup gatk GenomicsDBImport --batch-size 50 --sample-name-map ../info/all.map --genomicsdb-workspace-path chr08.database -L SL3.0ch08 >chr8.GDB.oe &
# nohup gatk GenomicsDBImport --batch-size 50 --sample-name-map ../info/all.map --genomicsdb-workspace-path chr09.database -L SL3.0ch09 >chr9.GDB.oe &
# nohup gatk GenomicsDBImport --batch-size 50 --sample-name-map ../info/all.map --genomicsdb-workspace-path chr10.database -L SL3.0ch10 >chr10.GDB.oe &
# nohup gatk GenomicsDBImport --batch-size 50 --sample-name-map ../info/all.map --genomicsdb-workspace-path chr11.database -L SL3.0ch11 >chr11.GDB.oe &
# nohup gatk GenomicsDBImport --batch-size 50 --sample-name-map ../info/all.map --genomicsdb-workspace-path chr12.database -L SL3.0ch12 >chr12.GDB.oe &

# nohup gatk GenotypeGVCFs -R ~/tomato/RAD2Core/genome/ref.normalized.fa -V gendb://chr01.database -O chr01.vcf.gz > chr1.GT.oe &
# nohup gatk GenotypeGVCFs -R ~/tomato/RAD2Core/genome/ref.normalized.fa -V gendb://chr02.database -O chr02.vcf.gz > chr2.GT.oe &
# nohup gatk GenotypeGVCFs -R ~/tomato/RAD2Core/genome/ref.normalized.fa -V gendb://chr03.database -O chr03.vcf.gz > chr3.GT.oe &
# nohup gatk GenotypeGVCFs -R ~/tomato/RAD2Core/genome/ref.normalized.fa -V gendb://chr04.database -O chr04.vcf.gz > chr4.GT.oe &
# nohup gatk GenotypeGVCFs -R ~/tomato/RAD2Core/genome/ref.normalized.fa -V gendb://chr05.database -O chr05.vcf.gz > chr5.GT.oe &
# nohup gatk GenotypeGVCFs -R ~/tomato/RAD2Core/genome/ref.normalized.fa -V gendb://chr06.database -O chr06.vcf.gz > chr6.GT.oe &
# nohup gatk GenotypeGVCFs -R ~/tomato/RAD2Core/genome/ref.normalized.fa -V gendb://chr07.database -O chr07.vcf.gz > chr7.GT.oe &
# nohup gatk GenotypeGVCFs -R ~/tomato/RAD2Core/genome/ref.normalized.fa -V gendb://chr08.database -O chr08.vcf.gz > chr8.GT.oe &
# nohup gatk GenotypeGVCFs -R ~/tomato/RAD2Core/genome/ref.normalized.fa -V gendb://chr09.database -O chr09.vcf.gz > chr9.GT.oe &
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
-O 2nd.combined.vcf.gz

## Index vcf
tabix -p vcf 2nd.combined.vcf.gz

## Filter variants
echo Vilter variants...
gatk VariantFiltration \
-R ~/tomato/RAD2Core/genome/ref.normalized.fa \
-V 2nd.combined.vcf.gz \
-filter "QUAL < 1000.0 || QD < 20.0 || FS > 60.0 || SOR > 3.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || AN < 768" \
--filter-name "my_snp_filter" \
-O 2nd.combined.filtered.vcf.gz

## Select SNP variants
echo Select variants...
gatk SelectVariants \
-R ~/tomato/RAD2Core/genome/ref.normalized.fa \
-V 2nd.combined.filtered.vcf.gz \
-select-type SNP \
--restrict-alleles-to BIALLELIC \
--exclude-filtered \
-O 2nd.combined.filtered.snp.vcf.gz

## Select INDEL variants
gatk SelectVariants -R ~/tomato/RAD2Core/genome/ref.normalized.fa -V 2nd.combined.filtered.vcf.gz -select-type INDEL --exclude-filtered -O 2nd.combined.filtered.indel.vcf.gz