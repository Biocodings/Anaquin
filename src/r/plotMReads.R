#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

data <- read.csv('%3%/%4%', row.names=1, sep='\t')
data <- data[data$reads != 0,]
data <- MetaQuin(seqs=row.names(data), input=log2(data$input), measured=log2(data$reads))

plotRead(data)