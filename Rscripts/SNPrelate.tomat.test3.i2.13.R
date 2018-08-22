### Analysis population snp data by CoreArray tools gdsfmt and SNPrelate
### Last modified by Dang, 13/12/2017.
### Ref: http://corearray.sourceforge.net/tutorials/SNPRelate/
# Load the R packages: gdsfmt and SNPRelate
library(gdsfmt)
library(SNPRelate)
library(ggplot2)
# load vcf and transform to gds
vcf.fn <- "/home/dang/tomato/test3/stacks.ref/i2.i3.stack/i2.i3.samples.snps.vcf"
snpgdsVCF2GDS(vcf.fn, "/home/dang/tomato/test3/stacks.ref/i2.i3.stack/i2.i3.samples.snps.gds", method="biallelic.only")
# Open the GDS file
genofile <- snpgdsOpen("/home/dang/tomato/test3/stacks.ref/i2.i3.stack/i2.i3.samples.snps.gds")
# Summary
snpgdsSummary(genofile)
# Take out genotype data for the first 3 samples and the first 5 SNPs
(g <- read.gdsn(index.gdsn(genofile, "genotype"), start=c(1,1), count=c(5,3)))
# Try different LD thresholds for sensitivity analysis
# Get all selected snp id
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2, autosome.only=FALSE)
snpset.id <- unlist(snpset)

## Run PCA
pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=2, autosome.only=FALSE)
# variance proportion (%)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))
# Get population information
pop_code <- read.table("/home/dang/tomato/test3/info/i2.i3.pop.txt", header=FALSE) 
colnames(pop_code) <- c("sample.id", "pop")
# make a data.frame
# assume the order of sample IDs is as the same as population codes
tab <- data.frame(sample.id = pca$sample.id,
                  pop = pop_code$pop,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)
# Draw
plot(tab$EV2, tab$EV1, col=as.integer(tab$pop), xlab="eigenvector 2", ylab="eigenvector 1")
legend("bottomright", legend=levels(tab$pop), pch="o", col=1:nlevels(tab$pop))
# Try ggplot
p <- ggplot(tab, aes(x=EV2, y=EV1, col=pop)) + geom_point() + geom_text(aes(label=tab$sample.id), hjust=0, vjust=0) 
p
# Plot the principal component pairs for the first four PCs
lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
pairs(pca$eigenvect[,1:4], col=tab$pop, labels=lbls)

## Parallel coordinates plot for the top principal components
library(MASS)
datpop <- pop_code$pop
parcoord(pca$eigenvect[,1:16], col=datpop)

## To calculate the SNP correlations between eigenvactors and SNP genotypes
# Get chromosome index
chr <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))
#chr_list <- c("Aradu.A01","Aradu.A02","Aradu.A03","Aradu.A04","Aradu.A05",
#              "Aradu.A06","Aradu.A07","Aradu.A08","Aradu.A09","Aradu.A10",
#              "Araip.B01","Araip.B02","Araip.B03","Araip.B04","Araip.B05",
#              "Araip.B06","Araip.B07","Araip.B08","Araip.B09","Araip.B10")
CORR <- snpgdsPCACorr(pca, genofile, eig.which=1:4)
savepar <- par(mfrow=c(3,1), mai=c(0.3, 0.55, 0.1, 0.25))
for (i in 1:3)
{
  dat1<-data.frame(CORR=abs(CORR$snpcorr[i,]),chr=chr)
  #subset_dat1<-dat1[dat1$chr%in%chr_list,]
  #subset_dat1<-droplevels(subset_dat1)
  plot(dat1$CORR, ylim=c(0,1), xlab="", ylab=paste("PC", i),
       col=dat1$chr, pch="+")
}
par(savepar)

## Identity-By-State Analysis
ibs <- snpgdsIBS(genofile, num.thread=2, autosome.only=FALSE)
# individulas in the same population are clustered together
pop.idx <- order(pop_code$pop)
image(ibs$ibs[pop.idx, pop.idx], col=terrain.colors(16))
# To perform multidimensional scaling analysis on the n??n matrix of genome-wide IBS pairwise distances
loc <- cmdscale(1 - ibs$ibs, k = 2)
x <- loc[, 1]; y <- loc[, 2]
race <- as.factor(pop_code$pop)
plot(x, y, col=race, xlab = "", ylab = "",
     main = "Multidimensional Scaling Analysis (IBS)")
legend("topleft", legend=levels(race), text.col=1:nlevels(race))
# To perform cluster analysis on the n??nn??n matrix of genome-wide IBS pairwise distances, and determine the groups by a permutation score
# Determine groups of individuals automatically
ibs.hc <- snpgdsHCluster(snpgdsIBS(genofile, num.thread=2, autosome.only=FALSE ))
rv <- snpgdsCutTree(ibs.hc)
plot(rv$dendrogram, leaflab="none", main="HapMap Phase II")
table(rv$samp.group)
# Determine groups of individuals by population information
rv2 <- snpgdsCutTree(ibs.hc, samp.group=as.factor(pop_code$pop))
plot(rv2$dendrogram, leaflab="none", main="HapMap Phase II")
legend("topright", legend=levels(race), col=1:nlevels(race), pch=19, ncol=4)
# Close the GDS file
closefn.gds(genofile)
