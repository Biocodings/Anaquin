#
# Anaquin - Sequin statistical analysis
#
# This R script was generated at %1%.
#
#    %2%
#

library(plyr)
library(Anaquin)

data <- read.csv('%3%/%4%', sep='\t')

# Group insertions and deletions
data$Mutation <- revalue(data$Mutation, c(Insertion = 'Indel', Deletion = 'Indel'))

# Only the sequin region
data <- data[data$Name != '-',]

# Construct unique identifiers
data$Unique <- paste(paste(data$Name, data$Pos, sep='_'), data$Type, sep='_')

plotROC(data$Unique, data$%5%, data$Mutation, data$Label, title='ROC (ranked by %5%)', legTitle=%6%, refGroup=%7%)
