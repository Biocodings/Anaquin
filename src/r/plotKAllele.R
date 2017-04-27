#
# Anaquin - Sequin statistical analysis. Version 1.1.1.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

data <- read.csv('%3%/%4%', row.names=1, sep='\t')
data <- data[data$ExpRef > 0,]
data <- data[data$ObsRef > 0 & data$ObsVar > 0,]

# Expected allele frequency (x-axis)
input <- log2(data$ExpVar / (data$ExpRef + data$ExpVar))

# Measured allele frequency (y-axis)
measured <- log2(data$ObsVar / (data$ObsRef + data$ObsVar))

plotLinear(row.names(data), input, measured, title='Allele Frequency', xlab='Expected Allele Frequency (log2)', ylab='Observed Allele Frequency (log2)', showLOQ=TRUE)