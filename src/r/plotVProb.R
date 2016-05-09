#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

# Change this for your significance level
sign <- 0.10

data <- read.csv('%3%/%4%', sep='\t')

data$name <- paste(data$seq, data$pos, sep='_')
data$name <- paste(data$name, data$type, sep='_')

data[data$label=='FP',]$eFold  <- 1
data[data$label=='FP',]$eAllele <- 1

data <- data[data$type=='SNP',] # Can be changed to indels
data <- VarQuin(seqs=data$name, ratio=as.factor(data$eFold), pval=data$pval, measured=data$eAllele)

plotAlleleP(data)