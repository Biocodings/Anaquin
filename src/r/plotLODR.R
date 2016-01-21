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

# Change to your chooden FDR rate
chosenFDR <- 0.1

plotLODR(data, choseFDR=0.1, shouldTable=FALSE, lvl='%5%', shouldBand=FALSE, yBreaks=c(1e-300, 1e-200, 1e-100, 1e-10, 1.00), locBand='local')