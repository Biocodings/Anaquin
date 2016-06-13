#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

data <- read.csv('%3%/%4%', row.name=1, sep='\t')

title <- '%5%'
xlab  <- '%6%'
ylab  <- '%7%'

expected <- %8%
measured <- %9%

data <- Anaquin(seqs=row.names(data), expected=expected, measured=measured)

plotScatter(data, title=title, xlab=xlab, ylab=ylab, showLOQ=%10%)



#data <- TransQuin(seqs=row.names(data), input=log2(data$input), measured=log2(data[,2:ncol(data)]))

plotScatter(data, title=title, xlab=xlab, ylab=ylab, showLOQ=%10%)