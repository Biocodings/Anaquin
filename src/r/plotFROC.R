#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

data <- read.csv('%3%/%4%', sep='\t')
data <- FusQuin(seqs=paste(data$pos1, data$pos2), expected=1, label=data$label, measured=data$read)
 
plotROC(data, title='ROC for FusQuin')