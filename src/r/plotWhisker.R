#
# Anaquin - Sequin statistical analysis
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

d1 <- read.csv('%3%/%4%', sep='\t')
d2 <- read.csv('%3%/%5%', sep='\t')
d3 <- read.csv('%3%/%6%', sep='\t')

plotWhisker(d1, d2, d3)