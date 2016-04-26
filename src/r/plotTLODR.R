#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

data <- read.csv('%3%/%4%', row.name=1, sep='\t')
data <- data[!is.na(data$expected),]
data <- TransQuin(seqs=row.names(data), mean=data$mean, expected=data$expected, measured=data$measured, pval=data$pval)

# Choose your FDR rate
chosenFDR <- 0.1

plotLODR(data, chosenFDR=chosenFDR)