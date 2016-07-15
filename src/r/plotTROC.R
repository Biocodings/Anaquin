#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

# Load the reference sequins
data <- read.csv('%3%/%4%', row.names=1, sep='\t')

# Create Anaquin data set
data <- CreateDataForAnaquin(names=row.names(data), expected=data$ExpLFC, measured=data$ObsLFC, score=1-data$Pval, qval=data$Qval)

data$seqs <- TransDiff_(data)

# Change this for another title
title <- 'ROC Plot'

plotROC(data, title=title, refRats=0)