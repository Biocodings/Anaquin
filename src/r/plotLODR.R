#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

setwd('%3%')

# Read the count table
counts <- read.csv('%4%', row.name=1)

# Read the differential results
diffs <- read.csv('%5%', row.name=1)

# Create a TransQuin data set for Anaquin
data <- transQuin(seqs=row.names(diffs), counts=rowMeans(counts), pval=diffs$pval)

# Change to your chooden FDR rate
chosenFDR <- 0.1

plotLODR(data, chosenFDR, lvl='%6%')