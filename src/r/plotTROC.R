#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

data <- read.csv('%3%/%4%', row.names=1, sep='\t')
data <- data[!is.na(data$Expected),]
data <- Anaquin(seqs=row.names(data), expected=data$Expected, measured=data$Measured, score=1-data$Pval, qval=data$Qval)
data$seqs <- TransDiff_(data)

title <- 'ROC Plot'

plotROC(data, title=title, refRatio=0)
