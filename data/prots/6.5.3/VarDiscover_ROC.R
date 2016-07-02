#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at 01-07-2016 10:16:28.
#
#    anaquin VarDiscover -o 6.5.3 -m MVA011.v013.csv -rvcf AVA009.v032.vcf -rbed AVA017.v001.bed -ufiles variants.txt 
#

library(Anaquin)

# Load the called variants
data <- read.csv('/data/Anaquin/6.5.3/6.5.3/VarDiscover_queries.csv', sep='\t')

# Remove false-positves that are not called within the sequin regions
data <- data[data$ID != '-',]

# Define your signifance
sign <- 0.10

if (all(data$Pval == '-'))
{
    data <- data
} else
{
    # It's only really FP if it's p-value is smaller than our threshold
    data <- data[data$Label!='FP' | (data$Pval == '-' | data$Pval<=sign),]
}

# Change to 'Indel' or delete the line for all variants
data <- data[data$Type=='SNP',]

# Change this for another title
title <- 'ROC Plot'

# Change this for another legend title
legTitle <- 'Allele Freq.'

# How to rank the ROC points
score <- 1-data$Pval

data$name <- paste(data$ID, data$Pos, sep='_')
data$name <- paste(data$name, data$Type, sep='_')

r <- unique(sort(data$EFold))

# Create Anaquin data set
data <- Anaquin(seqs=data$name, expected=data$EFold, score=score, label=data$Label)

plotROC(data, title=title, legTitle=legTitle, refRats=r)
