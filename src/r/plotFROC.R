#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

data <- read.csv('%3%/%4%', row.names=1, sep='\t')
data$name <- paste(data$pos1, data$pos2)
data <- FusQuin(seqs=data$name, expected=1, label=data$label, measured=data$reads)
 
plotROC(data, title='ROC for FusQuin')