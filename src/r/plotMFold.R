#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

#
# This script generates a log-fold plot for expected log-fold against measured log-fold.
#
#    - x-axis: expected log-fold
#    - y-axis: measured log-fold
#

library(Anaquin)

data <- read.csv('%3%/%4%', row.names=1, sep='\t')
data <- data[!is.na(data$ELFold),]
data <- TransQuin(seqs=row.names(data), input=data$ELFold, measured=data$MLFold)

plotTLogFold(data)