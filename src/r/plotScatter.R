#
# Anaquin - Sequin statistical analysis. Version 1.1.1.
#
# This R script was generated at %1%.
#
#    %2%
#

library('Anaquin')

# Load reference sequins
data <- read.csv('%3%/%4%', row.names=1, sep='\t')

title <- '%5%'
xlab  <- '%6%'
ylab  <- '%7%'

# Expected input concentration (x-axis)
input <- %8%

# Measured expression (y-axis)
measured <- %9%

# Create Anaquin data for plotScatter
anaquin <- createAnaquinData(names=row.names(data), %10%=input, measured=measured%12%)

plotScatter(anaquin, title=title, xlab=xlab, ylab=ylab, showLOQ=%11%)