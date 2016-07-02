#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

#
# Generate ROC curves for RNA-Seq differential analysis.
#

data <- read.csv('%3%/%4%', row.names=1, sep='\t')

# Create Anaquin data set
data <- Anaquin(seqs=row.names(data), expected=data$Expected, measured=data$Measured, score=1-data$Pval, qval=data$Qval)

data$seqs <- TransDiff_(data)

# Change this for another title
title <- 'ROC Plot'

plotROC(data, title=title, refRats=0) 