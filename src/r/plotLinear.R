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

# Expected input concentration (x-axis)
input <- %8%

# Measured expression (y-axis)
measured <- %9%

plotLinear(row.names(data), input, measured, title=title, xlab=xlab, ylab=ylab, showLOQ=%11%)