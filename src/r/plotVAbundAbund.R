#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

#
# Create a plot for expected abundance vs measured abundance
#

library(Anaquin)

data <- read.csv('%3%/%4%', row.names=1, sep='\t')
data <- VarQuin(seqs=row.names(data), expect=data$EAbund, measured=data$MAbund)

plotVAbundAbund(data)