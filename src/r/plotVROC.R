#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

# Choose your significance level
sign <- 0.10

data <- read.csv('%3%/%4%', sep='\t')
data <- data[!is.na(data$pval) & data$pval <= sign,]

data$name <- paste(data$sequin, data$pos, sep='_')
data$name <- paste(data$name, data$type, sep='_')

data <- VarQuin(seqs=data$name, expected=log2(data$eFold), pval=data$pval, label=data$label)

plotROC(data)