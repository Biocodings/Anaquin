#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

data <- read.csv('%3%/%4%', row.name=1, sep='\t')
data <- TransQuin(seqs=row.names(data), input=log2(data$input), measured=log2(data[,2:ncol(data)]))

plotExpress(data)