#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

#
# This script generates a plot for expected expression against measured expression.
#
#    - x-axis: expected expression in attomol/ul
#    - y-axis: measured expression
#

library(Anaquin)

data <- read.csv('%3%/%4%', row.names=1, sep='\t')
data <- VarQuin(seqs=row.names(data), expected=data$expected, measured=data$measured)

plotExpress(data)