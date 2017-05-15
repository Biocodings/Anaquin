#
# Anaquin - Sequin statistical analysis. Version 1.1.1.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

data <- read.csv('%3%/%4%', row.names=1, sep='\t')

input <- data$Expected
measured <- log2(data$Observed)

plotLinear(row.names(data), input, measured, title='CNV Ladder', xlab='Expected copy number', ylab='Observed Abundance (log2)', showLOQ=FALSE, showLinear=TRUE)