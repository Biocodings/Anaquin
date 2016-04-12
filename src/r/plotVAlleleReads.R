#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

#
# Create a plot for expected allele frequency vs measured read counts
#

library(Anaquin)

data <- read.csv('%3%/%4%', row.names=1, sep='\t')
data <- VarQuin(seqs=row.names(data), expected=data$expected, rRead=data$rCount, vRead=data$vCount, type=data$type)

plotAlleleReads(data)