#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

data <- read.csv('%3%/%4%', sep='\t')

# Choose your significance level
sign <- 0.10

# Change to 'Indel' or delete the line for all variants
data <- data[data$Type=='SNP',]

# Change this for another title
title <- 'ROC Plot'

# Change this for another legend title
legTitle <- 'Allele Freq.'

# Filter by p-values
data <- data[!is.na(data$Pval) & data$Pval <= sign,]

# Create Anaquin dataset
data <- Anaquin(seqs=paste(data$Seq, data$Pos, sep='_'), pval=data$Pval, label=data$Label, type=data$Type)

plotROC(data, title=title, legTitle=legTitle)
