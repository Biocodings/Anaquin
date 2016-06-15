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

expected <- %8%
measured <- %9%

data <- Anaquin(seqs=row.names(data), expected=expected, measured=measured)

plotScatter(data, title=title, xlab=xlab, ylab=ylab, showIntercept=%10%)