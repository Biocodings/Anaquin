#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library(plyr)
library(Anaquin)

data <- read.csv('%3%/%4%', sep='\t')

# Remove false positives that are not called within the sequin regions
data <- data[data$ID != '-',]

# How to rank the ROC points
score <- %5%

# Construct unique identifiers for the variants
data$Unique <- paste(paste(data$ID, data$Pos, sep='_'), data$Type, sep='_')

# Calculate the allele frequency
data$AlleleF <- round(data$ExpRef / data$ExpVar)

# Required for the next step
if (nrow(data[is.nan(data$AlleleF),]))
{
    data[is.nan(data$AlleleF),]$AlleleF <- 2
}

# Give better names for the groups
data$AlleleF <- revalue(as.factor(data$AlleleF), c('0'='Homozygous', '1'='Heterozygous', '2'='FP'))

# Create Anaquin data for PlotROC
anaquin <- AnaquinData(analysis='PlotROC', seqs=data$Unique, ratio=data$AlleleF, score=score, label=data$Label)

plotROC(anaquin, title=title, legTitle=legTitle, refRats=%6%)