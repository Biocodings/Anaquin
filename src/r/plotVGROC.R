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

data$FreqGrp <- revalue(as.factor(data$ExpFreq), c('0.500000'='Homozygous', '1.000000'='Heterozygous', '-'='FP'))
plotROC(data$Unique, data$Depth, data$FreqGrp, data$Label, title='ROC (ranked by VCF Depth)', legTitle='Allele Frequency', refGroup='FP')
