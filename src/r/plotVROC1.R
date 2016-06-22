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

# Load all the called variants
data <- read.csv('%3%/%4%', sep='\t')

fp <- data[data$Label=='FP',]
fp <- fp[with(fp, order(Pval)),]

data[data$Label=='FP',]$EFold <- 0

# Define your signifance
sign <- 0.10

# It's only really FP if it's p-value is smaller than our threshold
data <- data[data$Label!='FP' | data$Pval<=sign,]

# Change to 'Indel' or delete the line for all variants
data <- data[data$Type=='SNP',]

# Change this for another title
title <- 'ROC Plot'

# Change this for another legend title
legTitle <- 'Allele Freq.'

data$name <- paste(data$ID, data$Pos, sep='_')
data$name <- paste(data$name, data$Type, sep='_')

data <- Anaquin(seqs=data$name, expected=data$EFold, score=1-data$Pval, label=data$Label)

plotROC(data, title=title, legTitle=legTitle, refRatio=0)
