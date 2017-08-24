#
# Anaquin - Sequin statistical analysis
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

data <- read.csv('%3%/%4%', row.names=1, sep='\t')

# Expected allele frequency (x-axis)
input <- log2(data$ExpFreq)

# Measured allele frequency (y-axis)
measured <- log2(data$ObsFreq)

plotLinear(row.names(data), input, measured, title='Allele Frequency', xlab='Expected Allele Frequency (log2)', ylab='Observed Allele Frequency (log2)', showLOQ=FALSE)