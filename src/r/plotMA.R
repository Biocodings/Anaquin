#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

# Read the differential results
data <- read.csv('%3%/%4%', row.name=1)

# Create a TransQuin data set for Anaquin
data <- TransQuin(seqs=row.names(data), baseMean=data$baseMean, log2FoldChange=data$lfc, lfcSE=data$lfcSE, pvalue=data$pval, expected.LFC=data$elfc)

plotMA(data, lvl='%5%', shouldEndo=TRUE)