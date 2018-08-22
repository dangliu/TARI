# Use ggplot2

library(ggplot2)


# Read input

table <- read.table('/home/dang/tomato/test3/stacks.ref/i2.stack/SL3.popSNP.100kb.txt')
colnames(table) <- c('chr', 'start', 'end', 'count')
head(table)


# plot

p <- ggplot(data=table, aes(x=end/1000000, y=count))
p <- p + geom_line()
p <- p + facet_grid(chr ~ .)
p <- p + theme(axis.line = element_line(colour = "black"), panel.background = element_blank()) 
p <- p + labs(x="Chromosome position (Mb)", y="SNPs (n)")
p
