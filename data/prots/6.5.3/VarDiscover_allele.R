#
# Anaquin - Sequin statistical analysis. Version 1.1.1.
#
# This R script was generated at 01-07-2016 10:16:28.
#
#    anaquin VarDiscover -o 6.5.3 -m MVA011.v013.csv -rvcf AVA009.v032.vcf -rbed AVA017.v001.bed -ufiles variants.txt 
#

library(Anaquin)
data <- read.csv('/data/Anaquin/6.5.3/6.5.3/VarDiscover_sequins.csv', row.names=1, sep='\t')

title <- 'Allele Frequency'
xlab  <- 'Expected allele frequency (log2)'
ylab  <- 'Measured allele frequency (log2)'

expected <- log2(data$EFreq)
measured <- log2(data$MFreq)

data <- Anaquin(seqs=row.names(data), expected=expected, measured=measured)

plotScatter(data, title=title, xlab=xlab, ylab=ylab, showLOQ=TRUE)
