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

data[data$label=='FP',]$ratio  <- 1
data[data$label=='FP',]$allele <- 1

data <- data[data$type=='SNP',] # Change it to indel if required
data <- VarQuin(seqs=data$name, ratio=as.factor(data$ratio), pval=data$pval, measured=data$allele)

plotAlleleP(data)