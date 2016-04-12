#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

data <- read.csv('%3%/%4%', sep='\t', row.name=1)
data <- TransQuin(seqs=row.names(data), expect=data$EAbund, A1=data[,2], A2=data[,3], A3=data[,4])

plotMultiAbundAbund(data)