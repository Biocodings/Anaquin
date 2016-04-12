#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

data <- read.csv('%3%/%4%', row.names=1, sep='\t')
data <- data[!is.na(data$Ratio),]

# Create a data set for Anaquin
data <- LadQuin(seqs=row.names(data), expected=data$Expected, measured=data$Observed)

plotScatter(data, title='Expected ladders vs measured ladders')
