#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library('Anaquin')

# Load called variants
data <- read.csv('%3%/%4%', sep='\t')

# Remove false-positves that are not called within the sequin regions
data <- data[data$ID != 'NA',]

# Uncomment the line to filter for SNPs (similar for 'Insertion' and 'Deletion')
#data <- data[data$Type=='SNP',]

title <- 'ROC Plot'
legTitle <- 'Allele Freq.'

# How to rank the ROC points
score <- %5%

data$name <- paste(data$ID, data$Pos, sep='_')
data$name <- paste(data$name, data$Type, sep='_')

data$ExpFold <- round(data$ExpRef / data$ExpVar)

# Give the false-positives a dummy ratio group
if (nrow(data[data$Label=='FP',]) > 0)
{
    data[data$Label=='FP',]$ExpFold <- 0
}

# Create Anaquin data for PlotROC
anaquin <- AnaquinData(analysis='PlotROC', seqs=data$name, ratio=data$ExpFold, score=score, label=data$Label)

plotROC(anaquin, title=title, legTitle=legTitle, refRats=%6%)