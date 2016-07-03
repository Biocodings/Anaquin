#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at 03-07-2016 04:22:34.
#
#    anaquin RnaFoldChange -o 5.6.3 -m MTR004.v013.csv -ufiles DESeq2.txt 
#

library(Anaquin)

#
# Generate ROC curves for RNA-Seq differential analysis.
#

data <- read.csv('/data/Anaquin/5.6.3/5.6.3/RnaFoldChange_sequins.csv', row.names=1, sep='\t')

# Create Anaquin data set
data <- Anaquin(seqs=row.names(data), expected=data$Expected, measured=data$Measured, score=1-data$Pval, qval=data$Qval)

data$seqs <- TransDiff_(data)

# Change this for another title
title <- 'ROC Plot'

plotROC(data, title=title, refRats=0) 
