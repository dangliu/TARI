# Use ggplot2

library(ggplot2)


# Read input

table <- read.table('/home/tari/tomato/RAD2Core/gatk/GVCF/missing_readsN.txt')
colnames(table) <- c('sample_ID', 'reads', 'missing_count', 'missing_rate')
head(table)


# plot

p <- ggplot(data=table, aes(x=reads/1000000, y=missing_rate))
p <- p + geom_point(alpha=0.5)
p <- p + geom_vline(xintercept=0.4, color="red")
p <- p + theme(axis.line = element_line(colour = "black"), panel.background = element_blank()) 
p <- p + labs(x="Number of reads after QC (n*10^6)", y="Missing rate (%)")
p

