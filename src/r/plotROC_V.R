#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

# Change this for your significance level
sign <- 0.10

data <- read.csv('%3%/%4%', sep='\t')
data <- data[data$PValue <= sign,]

data$name <- paste(data$Sequin, data$Position, sep='_')
data$name <- paste(data$name, data$Type, sep='_')

# Create a VarQuin data set for Anaquin
data <- VarQuin(seqs=data$name, pval=data$PValue, expected=data$EAlleleF, label=data$Label, type=data$Type)

plotROC.VarQuin(data)
