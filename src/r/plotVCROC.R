#
# Anaquin - Sequin statistical analysis
#
# This R script was generated at %1%.
#
#    %2%
#

#
# Construct ROC curves for somatic mutation. This is a template, feel free to modify it.
#

library(plyr)
library(Anaquin)

data <- read.csv('%3%/%4%', sep='\t')

data <- data[data$Name != '-',]
data <- data[data$Context == 'Cancer' | data$Context == '-',]

# How to rank the ROC points
score <- %5%

# Construct unique identifiers for the variants
data$Unique <- paste(paste(data$Name, data$Pos, sep='_'), data$Type, sep='_')

# Recode the allele frequency
data$ExpFreq <- revalue(as.factor(data$ExpFreq), c('-'='FP'))

#
# 1. ROC for all true-positives
#

plotROC(data$Unique, score, data$Label, data$Label, title='ROC Plot (ranked by Allele Frequency)', legTitle='Sequins', refGroup=%6%)

#
# 2. ROC by zygosity
#

plotROC(data$Unique, score, data$ExpFreq, data$Label, title='ROC Plot (ranked by Allele Frequency)', legTitle='z', refGroup=%6%)