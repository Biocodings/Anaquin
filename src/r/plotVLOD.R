#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

data <- read.csv('%3%/%4%', sep='\t')

# False-positives have no probabilites
data <- data[data$Label != 'FP',]

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

data$EFold <- round(data$ERef / data$EVar)
data <- Anaquin(seqs=paste(data$Seq, data$Pos, sep='_'), expected=as.factor(data$EFold), pval=data$Pval, measured=data$EFreq)

plotLOD(data, title=title, xlab=xlab, ylab=ylab, legTitle=legTitle)  