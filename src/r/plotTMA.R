#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

data <- read.csv('%3%/%4%', row.name=1)
data <- TransQuin(seqs=row.names(data), mean=data$mean, expected=data$expected, measured=data$measured, se=data$se, pval=data$pval)

plotMA(data, lvl='%5%', shouldEndo=TRUE)