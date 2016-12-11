#!/usr/bin/Rscript
library(ggplot2)
dmelcontiglenths <- read.csv("../../data/processed/dmel-all-chromosome-r6.13.lengths",sep = "\t")
#ggplot(data=dmelcontiglenths, aes(Length)) + geom_histogram() + scale_x_log10()
ggplot(data=dmelcontiglenths, aes(Length)) + geom_histogram(bins = 100) + scale_x_log10(breaks=c(0,40000000,1000)) + labs(title = "Drosophila melanogaster Contig length Distribution", x = "Length", y = "Count (log scale)")
ggsave("../../output/figures/dmellenghts_histogram.png")
quit()
