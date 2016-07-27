#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library('Anaquin')

# Load reference sequins
data <- read.csv('%3%/%4%', row.name=1, sep='\t')

# Remove undetected sequins
data <- data[!is.na(data$ObsLFC),]

# Choose your FDR rate
FDR <- 0.1

xlab  <- 'Average Counts'
ylab  <- 'P-value'
title <- 'LODR Curves'

# Measured abundance
measured <- data$Mean

# Expected log-fold change
ratio <- data$ExpLFC

# Measured p-value
pval <- data$Pval

# Measured q-value
qval <- data$Qval

# Create Anaquin data for plotLODR
anaquin <- createAnaquinData(names=row.names(data), measured=measured, ratio=ratio, pval=pval, qval=qval)

plotLODR(anaquin, xlab=xlab, ylab=ylab, title=title, FDR=FDR, legTitle='LFC')