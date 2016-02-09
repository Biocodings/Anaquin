#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

data <- read.csv('%3%/%4%', row.names=1, sep='\t')

# Filter only the SNP
snp <- data[data$Type=='SNP',]

# Filter only the indels
ind <- data[data$Type=='Indel',]

#
# 1. Plot for all variants
#

data <- VarQuin(seqs=row.names(data), expected=data$EAlleleF, measured=data$MAlleleF)
plotScatter(data, title='Expected vs measured allele frequency (all variants)', xname='Expected log2 allele frequency', yname='Measured log2 allele frequency')

#
# 2. Plot for SNPs
#

data <- VarQuin(seqs=row.names(snp), expected=snp$EAlleleF, measured=snp$MAlleleF)
plotScatter(data, title='Expected vs measured allele frequency (SNP)', xname='Expected log2 allele frequency', yname='Measured log2 allele frequency')

#
# 2. Plot for indels
#

data <- VarQuin(seqs=row.names(ind), expected=ind$EAlleleF, measured=ind$MAlleleF)
plotScatter(data, title='Expected vs measured allele frequency (Indel)', xname='Expected log2 allele frequency', yname='Measured log2 allele frequency')