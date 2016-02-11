#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

data <- read.csv('%3%/%4%', row.names=1, sep='\t')

# Create a data set for Anaquin
data <- VarQuin(seqs=row.names(data), expected=data$EAlleleF, rRead=data$RCount, vRead=data$VCount, type=data$Type)

# Plot for SNV
plotVCoverage(data, title='SNV Allele Frequency')