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

data[data$Label=='FP',]$ERatio   <- 1
data[data$Label=='FP',]$EAlleleF <- 1

snp <- data[data$Type=='SNP',]
ind <- data[data$Type=='Indel',]

# Plot for SNP
plotLODR.Plot(data.frame(pval=log10(snp$PValue), abund=log10(snp$EAlleleF), ratio=as.factor(snp$EAlleleF)), title='Expected allele frequency vs P-value (SNP)', xname='Expected allele frequency (log10)', yname='P-value (log10)')

# Plot for indels
plotLODR.Plot(data.frame(pval=log10(ind$PValue), abund=log10(ind$EAlleleF), ratio=as.factor(ind$EAlleleF)), title='Expected allele frequency vs P-value (Indel)', xname='Expected allele frequency (log10)', yname='P-value (log10)')
