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

# Change to your chooden FDR rate
chosenFDR <- 0.1

xBreaks <- c(1, 10, 100, 1000, 10000)
xLabels <- c('1e+00', '1e+01', '1e+02', '1e+03', '1e+04')
yBreaks <- c(1e-300, 1e-200, 1e-100, 1e-10, 1.00)

plotLODR(data, shouldBand=TRUE, locBand='global', choseFDR=0.1, shouldTable=FALSE, lvl='gene', yBreaks=yBreaks, xBreaks=xBreaks, xLabels=xLabels)