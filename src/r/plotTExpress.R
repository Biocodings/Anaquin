#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

data <- read.csv('%3%/%4%', row.names=1, sep='\t')
data <- TransQuin(seqs=row.names(data), expected=log2(data$input), measured=log2(data$measured))

plotExpress(data)