#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

# Load reference sequins
data <- read.csv('%3%/%4%', row.name=1, sep='\t')

# Remove undetected sequins
data <- data[!is.na(data$ObsLFC),]

# Create Anaquin data set
data <- CreateDataForAnaquin(names=row.names(data), measured=data$Mean, ratio=data$ExpLFC, pval=data$Pval)

# Choose your FDR rate
fdr <- 0.1

# Label for the x-axis
xlab <- 'Average Counts'

# Label for the y-axis
ylab <- 'P-value'

# Title of the plot
title <- 'LODR Curves'

plotLODR(data, xlab=xlab, ylab=ylab, title=title, fdr=fdr)