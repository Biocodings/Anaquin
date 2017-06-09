#
# Anaquin - Sequin statistical analysis
#
# This R script was generated at %1%.
#
#    %2%
#

#
# Construct ROC curves for WGS. This is a template, feel free to modify it.
#

library(plyr)
library(Anaquin)

data <- read.csv('%3%/%4%', sep='\t')

# Only the sequin region (remove it for all variants)
data <- data[data$Name != '-',]

# How to rank the ROC points
score <- %5%

# Construct unique identifiers for the variants
data$Unique <- paste(paste(data$Name, data$Pos, sep='_'), data$Type, sep='_')

#
# 1. ROC for all sequin TP and FP 
#

plotROC(data$Unique, score, data$Label, data$Label, title='ROC (ranked by VCF Depth)', legTitle='Sequins', refGroup=%6%)

#
# 2. ROC by genotype
#

data$FreqGrp <- revalue(as.factor(data$ExpFreq), c('0.500000'='Homozygous', '1.000000'='Heterozygous', '-'='FP'))
plotROC(data$Unique, score, data$FreqGrp, data$Label, title='ROC (ranked by VCF Depth)', legTitle='Genotype', refGroup=%6%)

#
# 3. ROC by context
#

plotROC(data$Unique, score, data$Context, data$Label, title='ROC Plot (ranked by VCF Depth)', refGroup='-')