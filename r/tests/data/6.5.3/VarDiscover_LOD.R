#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at 01-07-2016 10:16:28.
#
#    anaquin VarDiscover -o 6.5.3 -m MVA011.v013.csv -rvcf AVA009.v032.vcf -rbed AVA017.v001.bed -ufiles variants.txt 
#

library(Anaquin)

data <- read.csv('/data/Anaquin/6.5.3/6.5.3/VarDiscover_queries.csv', sep='\t')

# False-positives have no probabilites
data <- data[data$Label != 'FP',]

# Change to 'Indel' or delete the line for all variants
data <- data[data$Type=='SNP',]

# Change this for another title
title <- 'Limit of Detection'

# Change this for another legend title
legTitle <- 'Allele Freq.'

# Change this for the x-axis label
xlab='Expected allele frequency (log10)'

# Change this for the y-axis label
ylab='P-value (log10)'

data <- Anaquin(seqs=paste(data$Seq, data$Pos, sep='_'), expected=as.factor(data$EFold), pval=data$Pval, measured=data$EAllele)

plotLOD(data, title=title, xlab=xlab, ylab=ylab, legTitle=legTitle)  
