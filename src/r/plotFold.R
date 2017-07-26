#
# Anaquin - Sequin statistical analysis
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

data <- read.csv('%3%/%4%', row.names=1, sep='\t')

title <- '%5%'
xlab  <- '%6%'
ylab  <- '%7%'

# Expected log-fold (x-axis)
expected <- %8%

# Measured log-fold (y-axis)
measured <- %9%%11%

plotLinear(row.names(data), expected, measured, title=title, xlab=xlab, ylab=ylab, showAxis=%10%, showLOQ=FALSE)
