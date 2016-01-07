#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

# Create a data set for Anaquin
data <- aqdata(seqs     = c(%3%),
               expected = c(%4%),
               measured = c(%5%))

plotScatter(data, xname = '%6%', yname = '%7%')