#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

# Change this for your significance level
signLevel <- 0.10

# Read the data file for true positives
tps <- read.csv('%3%/%4%.csv', sep='\t')

# Read the data file for false positives
fps <- read.csv('%3%/%5%.csv', sep='\t')

# Filter only the significant FPs
fps <- fps[fps$PValue < signLevel,]

tps$label <- 'TP'
fps$label <- 'FP'

data <- rbind(tps, fps)

data$name <- paste(data$Sequin, data$Position, sep='_')
data$name <- paste(data$name, data$Type, sep='_')

# Create a FusQuin data set for Anaquin
data <- FusQuin(seqs=data$name, pval=data$PValue, rRead=data$RefRead, vRead=data$VarRead, eAFreq=data$EAlleleF, label=data$label, type=data$Type)

plotROC.VarQuin(data)
