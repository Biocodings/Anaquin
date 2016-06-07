#
# Anaquin - Sequin statistical analysis. Version 1.1.1.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

title <- %7%
xlab  <- %8%
ylab  <- %9%

expected <- log2(data$%5%)
measured <- log2(data$%6%)

data <- read.csv('%3%/%4%', row.names=1, sep='\t')
data <- Anaquin(seqs=row.names(data), expected=expected, measured=measured)

plotScatter(data, title=title, xlab=xlab, ylab=ylab)