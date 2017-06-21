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

# Only the sequin region (remove it for all variants)
data <- data[data$Name != '-',]

# Construct unique identifiers for the variants
data$Unique <- paste(paste(data$Name, data$Pos, sep='_'), data$Type, sep='_')

#
# ROC for all sequin TP and FP 
#

plotROC(data$Unique, score, data$Label, data$Label, title='ROC (ranked by VCF Depth)', legTitle='Sequins', refGroup=%6%)

#
# ROC by genotype
#

data$FreqGrp <- as.factor(data$ExpFreq)
plotROC(data$Unique, %5%, data$FreqGrp, data$Label, title='ROC (ranked by VCF Depth)', legTitle='Allele Frequency', refGroup=%6%)