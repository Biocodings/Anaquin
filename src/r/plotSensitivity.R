#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

# Load the assembled sequins
data <- read.csv('%3%/%4%', row.names=1, sep='\t')

# Specify the sensitivity threshold
threshold <- 0.70

title <- '%5%'
xlab  <- '%6%'
ylab  <- '%7%'

# Create Anaquin data set
data <- Anaquin(seqs=row.names(data), expected=%8%, measured=%9%)

plotSensitivity(data, title=title, xlab=xlab, ylab=ylab, threshold=threshold, showLOA=%10%)