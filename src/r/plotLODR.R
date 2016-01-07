#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

# Create a data set for Anaquin
data <- aqdata(seqs   = c(%3%),
               counts = c(%4%),
               pval   = c(%5%),
               ratio  = c(%6%))

# Change to your chooden FDR rate
chosenFDR <- 0.1

plotLODR(data, chosenFDR)