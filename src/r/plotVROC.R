#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

# Load called variants
data <- read.csv('%3%/%4%', sep='\t')

# Remove false-positves that are not called within the sequin regions
data <- data[data$ID != 'NA',]

# Define your signifance
sign <- 0.10

if (all(is.na(data$Pval)))
{
    data <- data
} else
{
    # It's only really FP if it's p-value is smaller than our threshold
    data <- data[data$Label!='FP' | (data$Pval == '-' | data$Pval<=sign),]
}

# Change to 'SNP' or 'Indel'
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

# Create Anaquin data for plotROC
anaquin <- createAnaquinData(names=data$name, input=data$ExpFold, score=score, label=data$Label)

plotROC(anaquin, title=title, legTitle=legTitle, refRats=0)