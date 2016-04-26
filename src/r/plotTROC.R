#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

data <- read.csv('%3%/%4%', row.names=1, sep='\t')
data <- data[!is.na(data$expected),]
data <- TransQuin(seqs=row.names(data), expected=data$expected, measured=data$measured, pval=data$pval, qval=data$qval)
 
plotROC(data, title='ROC for TransQuin differential')