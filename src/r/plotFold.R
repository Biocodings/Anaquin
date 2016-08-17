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

# Expected log-fold (x-axis)
expected <- %8%

# Measured log-fold (y-axis)
measured <- %9%

# Create Anaquin data for PlotLinear
anaquin <- AnaquinData(analysis='PlotLinear', seqs=row.names(data), input=expected, measured=measured%11%)

plotLinear(anaquin, title=title, xlab=xlab, ylab=ylab, showAxis=%10%, showLOQ=FALSE)
