#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

#
# Plot ROC curve for somatic variants
#

library(Anaquin)

# Load the input variant file
data <- read.csv('%3%/%4%', sep='\t')

# Change to 'Indel' or delete the line for all variants
data <- data[data$Type=='SNP',]

# Change this for another title
title <- 'ROC Plot'

# Change this for another legend title
legTitle <- 'Allele Freq.'

data$name <- paste(data$ID, data$Pos, sep='_')
data$name <- paste(data$name, data$Type, sep='_')

data <- Anaquin(seqs=data$name, expected=data$EFold, score=1-data$Pval, label=data$Label)

plotROC(data, title=title, legTitle=legTitle)