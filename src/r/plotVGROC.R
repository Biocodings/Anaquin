#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

#
# Construct ROC curve for WGS by genotypes (SNPS and Indels). This is a template, feel free to modify it.
#

library(plyr)
library(Anaquin)

data <- read.csv('%3%/%4%', sep='\t')

# How to rank the ROC points
score <- %5%

# Construct unique identifiers for the variants
data$Unique <- paste(paste(data$ID, data$Pos, sep='_'), data$Type, sep='_')

# Calculate the allele frequency
data$ExpFreq <- revalue(as.factor(data$ExpFreq), c('0.500000'='Homozygous', '1.000000'='Heterozygous', '-'='FP'))

# Create Anaquin data for PlotROC
anaquin <- AnaquinData(analysis='PlotROC', seqs=data$Unique, ratio=data$ExpFreq, score=score, label=data$Label)

plotROC(anaquin, title='ROC Plot (ranked by VCF Depth)', legTitle='Genotype', refRats=%6%)