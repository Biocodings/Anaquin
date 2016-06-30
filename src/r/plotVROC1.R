#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

# Load the called variants
data <- read.csv('%3%/%4%', sep='\t')

# Remove false-positves that are not called within the sequin regions
data <- data[data$ID != '-',]

# Give a reference ratio to the FP
data[data$Label=='FP',]$EFold <- 0

# Define your signifance
sign <- 0.10

# It's only really FP if it's p-value is smaller than our threshold
data <- data[data$Label!='FP' | is.nan(data$Pval) | data$Pval<=sign,]

# Change to 'Indel' or delete the line for all variants
data <- data[data$Type=='SNP',]

# Change this for another title
title <- 'ROC Plot'

# Change this for another legend title
legTitle <- 'Allele Freq.'

# How to rank the ROC points
score <- %5%

data$name <- paste(data$ID, data$Pos, sep='_')
data$name <- paste(data$name, data$Type, sep='_')

# Create Anaquin data set
data <- Anaquin(seqs=data$name, expected=data$EFold, score=score, label=data$Label)

plotROC(data, title=title, legTitle=legTitle, refRats=0)