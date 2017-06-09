#
# Anaquin - Sequin statistical analysis
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

data <- read.csv('%3%/%4%', row.names=2, sep='\t')

# Expected input concentration (x-axis)
input <- %8%

# Measured expression (y-axis)
measured <- %9%
# %11%
plotConjoint(data$Sequin, row.names(data), input, measured, title='%5%', xlab='%6%', ylab='%7%')