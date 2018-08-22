### Simulation by SimRAD
### Last modified by Dang, 16/01/2018.
### Ref: https://cran.r-project.org/web/packages/SimRAD/SimRAD.pdf
# Load the R packages: SimRAD and ggplot2
library(SimRAD)
library(ggplot2)

# Sequence input
# importing the Fasta file and sub selecting 25% of the contigst
rfsq <- ref.DNAseq("/home/dang/tomato/loci_predict/ref.fa", subselect.contigs = TRUE, prop.contigs = 0.25)

# length of the reference sequence:
width(rfsq)
# ratio for the cross-multiplication of the number of fragments and loci at the genomes scale:
genome.size <-828076960 # genome size: 828.1Mb
ratio <- genome.size/width(rfsq)
ratio
# computing GC content:
require(seqinr)
GC(s2c(rfsq))
# simulating random generated DNA sequence with characteristics equivalent to
#     the sub-selected reference genome for comparison purpose:
smsq <- sim.DNAseq(size=width(rfsq), GCfreq=GC(s2c(rfsq)))




### a double digestion (ddRAD)

# Restriction Enzyme 1
# Ref: https://www.neb.com/products/r0140-psti#Product%20Information
# Ref2: https://www.neb.com/products/r0525-msei#Product%20Information
#NsiI
cs_5p1 <- "ATGCA"
cs_3p1 <- "T"
# Restriction Enzyme 2
#MseI 
cs_5p2 <- "T"
cs_3p2 <- "TAA"

# Use ref to test
rfsq.dig <- insilico.digest(rfsq, cs_5p1, cs_3p1, cs_5p2, cs_3p2, verbose=TRUE)
rfsq.sel <- adapt.select(rfsq.dig, type="AB+BA", cs_5p1, cs_3p1, cs_5p2, cs_3p2)
# wide size selection (200-270):
wid.rfsq <- size.select(rfsq.sel,  min.size = 200, max.size = 270, graph=TRUE, verbose=TRUE)
# narrow size selection (210-260):
nar.rfsq <- size.select(rfsq.sel,  min.size = 210, max.size = 260, graph=TRUE, verbose=TRUE)
# the resulting fragment characteristics can be further examined:
boxplot(list(width(rfsq.sel), width(wid.rfsq), width(nar.rfsq)), names=c("All fragments",
                                                                         "Wide size selection", "Narrow size selection"), ylab="Locus size (bp)")

# Use sim to test
smsq.dig <- insilico.digest(smsq, cs_5p1, cs_3p1, cs_5p2, cs_3p2, verbose=TRUE)
smsq.sel <- adapt.select(smsq.dig, type="AB+BA", cs_5p1, cs_3p1, cs_5p2, cs_3p2)
# wide size selection (200-270):
wid.smsq <- size.select(smsq.sel,  min.size = 200, max.size = 270, graph=TRUE, verbose=TRUE)
# narrow size selection (210-260):
nar.smsq <- size.select(smsq.sel,  min.size = 210, max.size = 260, graph=TRUE, verbose=TRUE)
# the resulting fragment characteristics can be further examined:
boxplot(list(width(smsq.sel), width(wid.smsq), width(nar.smsq)), names=c("All fragments",
                                                                         "Wide size selection", "Narrow size selection"), ylab="Locus size (bp)")
