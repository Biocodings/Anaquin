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
data <- data[data$pval <= sign,]

data$name <- paste(data$sequin, data$Position, sep='_')
data$name <- paste(data$name, data$Type, sep='_')

data <- VarQuin(seqs=data$name, expected=data$expected, pval=data$pval, label=data$label, type=data$type)

plotROC(data)