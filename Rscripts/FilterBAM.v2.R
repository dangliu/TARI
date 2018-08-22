#!/usr/bin/Rscript

# check package
library(Rsamtools)
library(reshape2)

args<-commandArgs(TRUE)

qlist <- function(bam.in, qname.list= "IncludeQueryName.txt"){
	filterParm  <- ScanBamParam(what = c("qname","rname","flag", "pos"), tag = c("XS", "YT", "DT"))
	ftag= scanBam(bam.in, param= filterParm) # sort by name and compress
	if (is.null(ftag[[1]]$tag$DT)){ftag[[1]]$tag$DT <- rep(NA, length(ftag[[1]]$qname))} # Modify for no duplicates
	x <- DataFrame(ftag[[1]])
	isRead1 <- bamFlagAsBitMatrix(x$flag)[,"isFirstMateRead"]==1L
	isRead2 <- !isRead1
	x.pool= merge(x[isRead1, ], x[isRead2, ], by= "qname", all= TRUE, sort= FALSE, suffixes= c(".R1", ".R2"))
	key= (
		is.na(x.pool$"tag.XS.R1") & is.na(x.pool$"tag.XS.R2") & # unique map
		x.pool$"tag.YT.R1"== "CP" & # concordant map
		is.na(x.pool$"tag.DT.R1") & # remove optical duplicate DT.tag== "SQ"
		as.character(x.pool$"rname.R1") %in% c("SL3.0ch01", "SL3.0ch02", "SL3.0ch03", "SL3.0ch04", "SL3.0ch05", "SL3.0ch06", "SL3.0ch07", "SL3.0ch08", "SL3.0ch09", "SL3.0ch10", "SL3.0ch11", "SL3.0ch12") # & # only target chromosome
		)
	qnames= x.pool$"qname"[key]
	write.table(qnames, file= qname.list, row.names= FALSE, col.names= FALSE, quote= FALSE)
}


### build filter query name list
qlist(bam.in=args[1], qname.list=args[2])