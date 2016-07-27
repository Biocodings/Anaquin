#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library('Anaquin')

# Load reference sequins
data <- read.csv('%3%/%4%', row.names=1, sep='\t')

# Specify sensitivity threshold
threshold <- 0.70

title <- '%5%'
xlab  <- '%6%'
ylab  <- '%7%'

# Expected input concentration (x-axis)
input <- %8%

# Measured sensitivity (y-axis)
sn <- %9%

# Create Anaquin data for plotSensitivity
anaquin <- createAnaquinData(names=row.names(data), input=input, sensitivity=sn)

plotSensitivity(anaquin, title=title, xlab=xlab, ylab=ylab, threshold=threshold, showLOA=%10%)