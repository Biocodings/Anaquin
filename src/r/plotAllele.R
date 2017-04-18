#
# Anaquin - Sequin statistical analysis. Version 1.1.1.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

data <- read.csv('%3%/%4%', row.names=1, sep='\t')
data <- data[data$Label == 'TP' & data$Group == 'Cosmic',]

# Convert from factor representation
data$ObsFreq <- as.numeric(as.character(data$ObsFreq))

title <- '%5%'
xlab  <- '%6%'
ylab  <- '%7%'

# Expected input concentration (x-axis)
input <- %8%

# Measured expression (y-axis)
measured <- %9%

plotLinear(row.names(data), input, measured, title=title, xlab=xlab, ylab=ylab, showLOQ=%11%)