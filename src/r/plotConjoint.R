#
# Anaquin - Sequin statistical analysis
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

data <- read.csv('%3%/%4%', row.names=2, sep='\t')

# Expected copy number (x-axis)
x <- %8%

# Measured expression (y-axis)
y <- %9%

plotConjoint(data$Name, row.names(data), x, y, title='%5%', xlab='%6%', ylab='%7%')