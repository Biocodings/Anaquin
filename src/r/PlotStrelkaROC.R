#
# Anaquin - Sequin statistical analysis
#
# This R script was generated at %1%.
#
#    %2%
#

library(plyr)
library(Anaquin)

data <- read.csv('%3%/%4%', sep='\t')

data <- data[data$Name != '-',]
data <- data[data$Context == 'Cancer' | data$Context == '-',]

data$Score <- suppressWarnings(as.factor(ifelse(data$QSI == '-', as.numeric(as.character(data$QSS)), as.numeric(as.character(data$QSI)))))

# Unique identifiers for the variants
data$Unique <- paste(paste(data$Name, data$Pos, sep='_'), data$Mutation, sep='_')

plotROC(data$Unique, data$Score, data$ExpFreq, data$Label, title='ROC Plot (ranked by QSI/QSS)', legTitle='Allele Frequency', refGroup='-')
