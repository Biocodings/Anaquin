#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

data <- read.csv('%3%/%4%', row.names=1, sep='\t')
data <- data[!is.na(data$Ratio),]
data <- LadQuin(seqs=row.names(data), input=data$input, measured=data$measured, ratio=data$ratio)

plotNorm(data)