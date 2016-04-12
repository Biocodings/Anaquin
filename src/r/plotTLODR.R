#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

#
# This script generates a LODR (limit-of-detection-ratio) plot.
#

library(Anaquin)

data <- read.csv('%3%/%4%', row.name=1, sep='\t')
data <- data[!is.na(data$expected),]
data <- TransQuin(seqs=row.names(data), mean=data$mean, expected=data$expected, measured=data$measured, se=data$se, pval=data$pval)

# Choose your FDR rate
chosenFDR <- 0.1

plotLODR(data, chosenFDR=0.1)