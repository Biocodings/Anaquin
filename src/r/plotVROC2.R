#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

data <- read.csv('%3%/%4%', sep='\t')
data <- data[data$type=='SNP',]

data$name <- paste(data$sequin, data$pos, sep='_')
data$name <- paste(data$name, data$type, sep='_')

m <- max(data$ref + data$var)
data <- VarQuin(seqs=data$name, expected=1, pval=m-(data$ref + data$var), label=data$label)

plotROC(data)