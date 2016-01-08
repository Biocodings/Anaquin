#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

# Create a data set for Anaquin
data <- aqData(seqs = c(%3%),
               pval = c(%4%))

plotROC(data)