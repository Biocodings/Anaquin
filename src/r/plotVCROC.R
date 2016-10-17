#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

data <- read.csv('%3%/%4%', sep='\t')

# Remove false positives that are not called within the sequin regions
data <- data[data$ID != 'NA',]

# How to rank the ROC points
score <- %5%

# FP don't have a sequin ID. We'll need to construct unique identifiers.
data$unique <- as.factor(seq(1:length(data$ID)))#paste(paste(data$ID, data$Pos, sep='_'), data$Type, sep='_')

data$AlleleF <- round(data$ExpRef / data$ExpVar)

# Give dummy allele frequency for FP
if (nrow(data[is.nan(data$AlleleF),]))
{
    data[is.nan(data$AlleleF),]$AlleleF <- -1
}

# Create Anaquin data for PlotROC
anaquin <- AnaquinData(analysis='PlotROC', seqs=data$unique, ratio=data$AlleleF, score=score, label=data$Label)

plotROC(anaquin, title='ROC Plot', legTitle='Allele Freq.', refRats=%6%)