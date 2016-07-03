#
# Anaquin - Sequin statistical analysis. Version 1.1.1.
#
# This R script was generated at 03-07-2016 04:22:33.
#
#    anaquin RnaFoldChange -o 5.6.3 -m MTR004.v013.csv -ufiles DESeq2.txt 
#

library(Anaquin)
data <- read.csv('/data/Anaquin/5.6.3/5.6.3/RnaFoldChange_sequins.csv', row.names=1, sep='\t')

title <- 'Fold Change'
xlab  <- 'Expected fold change (log2)'
ylab  <- 'Measured fold change (log2)'

expected <- data$Expected
measured <- data$Measured

# Create Anaquin data set
data <- Anaquin(seqs=row.names(data), expected=expected, measured=measured)

plotFold(data, title=title, xlab=xlab, ylab=ylab, showIntercept=TRUE)

