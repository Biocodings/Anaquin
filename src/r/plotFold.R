#
# Anaquin - Sequin statistical analysis. Version 1.1.1.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

# Load referebce sequins
data <- read.csv('%3%/%4%', row.names=1, sep='\t')

title <- '%5%'
xlab  <- '%6%'
ylab  <- '%7%'

expected <- %8%
measured <- %9%

# Create Anaquin data set
data <- CreateDataForAnaquin(names=row.names(data), expected=expected, measured=measured)

plotScatter(data, title=title, xlab=xlab, ylab=ylab, showInter=%10%, showLOQ=FALSE)