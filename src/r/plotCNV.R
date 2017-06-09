#
# Anaquin - Sequin statistical analysis
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

data <- read.csv('%3%/%4%', row.names=1, sep='\t')

# Expected copy number on the x-axis
CNV <- data$Copy

# Measured abundance on the y-axis
measured <- data$After

plotLinear(row.names(data), CNV, measured, title='CNV Ladder', xlab='Expected copy number', ylab='Observed Abundance', showLOQ=FALSE, showLinear=TRUE)

data <- data[data$Copy == 2,]
plotLinear(row.names(data), log2(data$After), log2(data$Genome), title='Copy number 2n sequins and their genomic regions', xlab='Sequin coverage (log2)', ylab='Genomic coverage (log2)', showLOQ=FALSE, showLinear=TRUE)
