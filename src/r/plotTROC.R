#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

#
# This script generates an ROC curve for the TransQuin sequins
#

library(Anaquin)

data <- read.csv('%3%/%4%', row.names=1, sep='\t')
data <- data[!is.na(data$ELFold),]
data <- TransQuin(seqs=row.names(data), expect=data$ELFold, measured=data$MLFold, pval=data$PValue, qval=data$QValue,
                  ratio=abs(data$ELFold))

plotTROC(data)