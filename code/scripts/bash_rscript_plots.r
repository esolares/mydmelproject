#!/usr/bin/Rscript
library(ggplot2)

#generate CDF of  contig lengths
dmelcontiglenths <- read.csv("../../data/processed/dmel-all-chromosome-r6.13.lengths",sep = "\t")
#ggplot(data=dmelcontiglenths, aes(Length)) + geom_histogram() + scale_x_log10()
ggplot(data=dmelcontiglenths, aes(Length)) + geom_histogram(bins = 100) + scale_x_log10(breaks=c(0,40000000,1000)) + labs(title = "Drosophila melanogaster Contig length Distribution", x = "Length", y = "Count (log scale)")
ggsave("../../output/figures/dmel_lenghts_CDF.png")

#generate histogram of GC content
dmelgccontent <- read.csv("../../data/processed/dmel-all-chromosome-r6.13_gc_content.txt", header = FALSE, sep = "\t")
ggplot(data = dmelgccontent, aes(V2)) + geom_histogram() + labs(title = "Drosophila melanogaster GC Content per Contig Histogram", x = "GC Content %", y = "Count")
ggsave("../../output/figures/dmel_gc_histogram.png")

#generate genelength histogram
dmelgenelengths <- read.csv("../../data/processed/dmel-all-r6.13_genelengths.txt", header = FALSE, sep = "\t")
ggplot(data = dmelgenelengths, aes(V2)) + geom_histogram() + labs(title = "Drosophila melanogaster Gene Lengths Histogram", x = "Gene Lengths", y = "Count") + scale_x_log10()
ggsave("../../output/figures/dmel_genelengths_histogram.png")

#generate exon length histogram
dmelgenelengths <- read.csv("../../data/processed/dmel-all-r6.13_exonlengths.txt", header = FALSE, sep = "\t")
ggplot(data = dmelgenelengths, aes(V2)) + geom_histogram() + labs(title = "Drosophila melanogaster Exon Lengths Histogram", x = "Exon Lengths", y = "Count") + scale_x_log10()
ggsave("../../output/figures/dmel_exonlengths_histogram.png")

quit()
