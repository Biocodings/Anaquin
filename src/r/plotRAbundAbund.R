#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

# Read table of FPKM
data <- read.csv('%3%/%4%', row.name=1)

data <- data[data$A1!=0,] # TODO: Fix this
data <- data[data$A2!=0,]
data <- data[data$A3!=0,]

# Create a TransQuin data set for Anaquin
data <- TransQuin(seqs = row.names(data), expected=data$expected, A1=data$A1, A2=data$A2, A3=data$A3)

plotPoolScatter(data, xname = 'Expected concentration (log2 attomol/ul)', yname = 'Measured coverage (log2 FPKM)')