#
# Anaquin - Sequin statistical analysis. Version 1.1.1.
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

# Expected input concentration (the x-axis)
input <- %8%

# Measured expression (the y-axis)
measured <- %9%

# Create Anaquin data for plotScatter
data <- CreateDataForAnaquin(seqs=row.names(data), expected=input, measured=measured)

plotScatter(data, title=title, xlab=xlab, ylab=ylab, showLOQ=%10%)