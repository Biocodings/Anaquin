#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at 03-07-2016 04:22:34.
#
#    anaquin RnaFoldChange -o 5.6.3 -m MTR004.v013.csv -ufiles DESeq2.txt 
#

library(Anaquin)

data <- read.csv('/data/Anaquin/5.6.3/5.6.3/RnaFoldChange_sequins.csv', row.name=1, sep='\t')
data <- data[!is.na(data$Expected),]
data <- Anaquin(seqs=row.names(data), mean=data$Mean, input=data$Expected, measured=data$Measured, pval=data$Pval)

# Choose your FDR rate
chosenFDR <- 0.1

plotLODR(data, chosenFDR=chosenFDR)  
