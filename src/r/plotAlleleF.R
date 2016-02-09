#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

data <- read.csv('%3%/%4%', row.names=1, sep='\t')

# Create a TransQuin data set for Anaquin
data <- VarQuin(seqs=row.names(data), expected=data$EAlleleF, measured=data$MAlleleF)

plotScatter(data, xname='Expected log2 allele frequency', yname='Measured log2 allele frequency')