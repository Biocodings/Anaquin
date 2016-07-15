#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

# Load reference sequins
data <- read.csv('%3%/%4%', row.name=1, sep='\t')

# Remove undetected sequins
data <- data[!is.na(data$ObsLFC),]

# Create Anaquin data set
data <- CreateDataForAnaquin(names=row.names(data), mean=data$Mean, input=data$ExpLFC, pval=data$Pval)

# Choose your FDR rate
chosenFDR <- 0.1

plotLODR(data, chosenFDR=chosenFDR) 