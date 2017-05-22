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
measured <- data$Observed

plotLinear(row.names(data), CNV, measured, title='CNV Ladder', xlab='Expected copy number', ylab='Observed Abundance', showLOQ=FALSE, showLinear=TRUE)