#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

data <- read.csv('%3%/%4%', row.names=1)
data <- TransQuin(seqs=row.names(data), B1=data$B1, B2=data$B2, B3=data$B3)

plotMajor(data)
