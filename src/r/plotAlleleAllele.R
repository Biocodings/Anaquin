#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

#
# Create a plot for expected allele fraction vs measured allele fraction
#

library(Anaquin)

data <- read.csv('%3%/%4%', row.names=1, sep='\t')
data <- VarQuin(seqs=row.names(data), expect=data$EAlleleF, measured=data$MAlleleF)

plotAlleleAllele(data)
