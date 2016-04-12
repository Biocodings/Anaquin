#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

#
# This script generates a log-fold plot for expected log-fold against measured log-fold.
#
#    - x-axis: expected log-fold
#    - y-axis: measured log-fold
#

library(Anaquin)

data <- read.csv('%3%/%4%', row.name=1, sep='\t')
data <- data[!is.na(data$ELFold),]
data <- TransQuin(seqs=row.names(data), baseMean=data$baseMean, log2FoldChange=data$lfc, lfcSE=data$lfcSE, pvalue=data$pval, expected.LFC=data$elfc)

# Choose your FDR rate
chosenFDR <- 0.1

plotLODR(data, chosenFDR=0.1)