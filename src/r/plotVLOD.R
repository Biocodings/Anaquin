#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

data <- read.csv('%3%/%4%', sep='\t')

# Don't lose those zero p-values
data[data$Pval==0,]$Pval <- min(data[data$Pval>0,]$Pval)

# False-positives have no probability
data <- data[data$Label!='FP',]

# Change to 'Indel' or delete the line for all variants
data <- data[data$Type=='SNP',]

xlab     <- 'Expected allele frequency (log10)'
ylab     <- 'P-value (log10)'
title    <- 'LOD Plot'
legTitle <- 'Allele Freq'

data$EFold <- round(data$ExpRef / data$ExpVar)

xBreaks <- c(1e-3, 1e-2, 1e-1)                    
xLabels <- c('-3', '-2', '-1')
yBreaks <- c(1, 1e-100, 1e-200, 1e-300)
yLabels <- c(1, -100, -200, -300)

# Create Anaquin data for PlotLODR
anaquin <- AnaquinData(analysis='PlotLODR', seqs=paste(data$ID, data$Position, sep='_'), ratio=data$EFold, pval=data$Pval, measured=data$ExpFreq)

plotLODR(anaquin, shouldFit=FALSE, xBreaks = xBreaks, xLabels=xLabels, yBreaks=yBreaks, yLabels=yLabels, title=title, xlab=xlab, ylab=ylab, legTitle=legTitle)