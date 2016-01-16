#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

# Create a TransQuin data set for Anaquin
data <- transQuin(seqs = c(%3%),
                  pval = c(%4%))

# Classify whether they're TP or FP
data <- transClassify(data)

plotROC(data)