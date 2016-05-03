#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

data <- read.csv('%3%/%4%', row.names=1, sep='\t')
data <- FusQuin(seqs=row.names(data), expected=log2(data$expected), measured=data$measured)
 
plotSigmoid(data, title='Sensitivity Plot', xlab='Input concentration (log2)', ylab='Sensitivity')