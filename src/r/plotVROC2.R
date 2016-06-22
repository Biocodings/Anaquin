#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

#
# Plot ROC curve for germline variants
#

library(Anaquin)

data <- read.csv('%3%/%4%', sep='\t')

# Change this for another title
title <- 'ROC Plot'

# Change this to indels or all variants
data <- data[data$Type=='SNP',]

# Change this to not use the depth coverage for ranking curve points
score <- data$Depth

data$name <- paste(data$ID, data$Pos, sep='_')
data$name <- paste(data$name, data$Type, sep='_')

data <- Anaquin(seqs=data$name, expected=1, score=score, label=data$Label)

plotROC(data, title=title)  
