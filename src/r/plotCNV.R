#
# Anaquin - Sequin statistical analysis. Version 1.1.1.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

data <- read.csv('%3%/%4%', row.names=1, sep='\t')

# Expected copy number on the x-axis
CNV <- data$CNV

# Measured abundance on the y-axis
measured <- log2(data$Observed)

plotLinear(row.names(data), input, measured, title='CNV Ladder', xlab='Expected copy number', ylab='Observed Abundance (log2)', showLOQ=FALSE, showLinear=TRUE)