#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

data <- read.csv('%3%/%4%', row.names=1, sep='\t')
data <- FusQuin(seqs=row.names(data), expected=log2(data$expected), measured=log2(data$measured))
 
plotFold(data, title='Expected log-ratio vs Measured log-ratio', xlab='Expected ratio (log2)', ylab='Measured ratio (log2)')