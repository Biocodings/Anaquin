#
# Anaquin - Sequin statistical analysis
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

data <- read.csv('%3%/%4%', row.names=1, sep='\t')
data <- data[data$Context == 'Cancer' & data$Label == 'TP',]

# Convert from factor representation
data$ObsFreq <- as.numeric(as.character(data$ObsFreq_Tumor))

title <- '%5%'
xlab  <- '%6%'
ylab  <- '%7%'

# Expected allele frequency (x-axis)
input <- %8%

# Measured allele frequency (y-axis)
measured <- %9%

plotLinear(row.names(data), input, measured, title=title, xlab=xlab, ylab=ylab, showLOQ=%11%)