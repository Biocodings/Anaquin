#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

#
# Create a plot for expected allele fraction vs measured allele fraction
#

library(Anaquin)

data <- read.csv('%3%/%4%', row.names=1, sep='\t')

# Create a data set for all variants
allVars <- VarQuin(seqs=row.names(data), expect=data$EAlleleF, measured=data$MAlleleF)

plotAlleleAllele(allVars)

#
# Uncomment the code if you want to plot for only SNPs.
#

# snp <- data[data$Type=='SNP',]
# allSNPs <- VarQuin(seqs=row.names(snp), expected=snp$EAlleleF, measured=snp$MAlleleF)
# plotAlleleAllele(allSNPs)

#
# Uncomment the code if you want to plot for only indels.
#

# ind <- data[data$Type=='Indel',]
# allInds <- VarQuin(seqs=row.names(ind), expected=ind$EAlleleF, measured=ind$MAlleleF)
# plotAlleleAllele(allInds	)
