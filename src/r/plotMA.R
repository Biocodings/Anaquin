#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

# Read the count table
file <- read.csv('%3%', row.name=1)

# Create a TransQuin data set for Anaquin
data <- transQuin(seqs = row.names(data), A1=data$A1, A2=data$A2, A3=data$A3, B1=data$B1, B2=data$B2, B3=data$B3)

plotMA(data, lvl='%4%', shouldEndo=TRUE)
