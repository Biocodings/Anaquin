#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

# Load reference sequins
data <- read.csv('%3%/%4%', row.names=1, sep='\t')

# Create Anaquin data set
data <- createAnaquinData(names=row.names(data), input=data$ExpLFC, measured=data$ObsLFC, score=1-data$Pval, qval=data$Pval)

data$seqs <- TransDiff_(data)

# Title of the ROC plot
title <- 'ROC Plot'

plotROC(data, title=title, refRats=0)
