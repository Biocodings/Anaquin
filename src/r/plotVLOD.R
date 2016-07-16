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

# Don't lose those zero p-values
data[data$Pval==0,]$Pval <- min(data[data$Pval>0,]$Pval)

# False-positives have no probability
data <- data[data$Label!='FP',]

# Change to 'Indel' or delete the line for all variants
data <- data[data$Type=='SNP',]

# Change this for another title
title <- 'LODR Plot'

# Change this for another legend title
legTitle <- 'Allele Freq'

# Change this for the x-axis label
xlab='Expected allele frequency (log10)'

# Change this for the y-axis label
ylab='P-value (log10)'

data$EFold <- round(data$ExpRef / data$ExpVar)

# Create Anaquin data set
data <- CreateDataForAnaquin(names=paste(data$ID, data$Position, sep='_'), ratio=data$EFold, pval=data$Pval, measured=data$ExpFreq)

xBreaks=c(1e-3, 1e-2, 1e-1)                    
xLabels=c('-3', '-2', '-1')
yBreaks=c(1, 1e-100, 1e-200, 1e-300)
yLabels=c(1, -100, -200, -300)

plotLODR(data, shouldFit=FALSE, xBreaks = xBreaks, xLabels=xLabels, yBreaks=yBreaks, yLabels=yLabels, title=title, xlab=xlab, ylab=ylab, legTitle=legTitle)
