#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

# Load reference sequins
data <- read.csv('%3%/%4%', sep='\t')

# False-positives have no probability
data <- data[data$Label!='FP',]

# Change to 'Indel' or delete the line for all variants
data <- data[data$Type=='SNP',]

# Change this for another title
title <- 'Limit of Detection'

# Change this for another legend title
legTitle <- 'Allele Freq.'

# Change this for the x-axis label
xlab='Expected allele frequency (log10)'

# Change this for the y-axis label
ylab='P-value (log10)'

data$EFold <- round(data$ExpRef / data$ExpVar)

# Create Anaquin data set
data <- CreateDataForAnaquin(names=paste(data$ID, data$Position, sep='_'), ratio=data$EFold, pval=data$Pval, measured=data$ExpFreq)

plotLOD(data, title=title, xlab=xlab, ylab=ylab, legTitle=legTitle)
