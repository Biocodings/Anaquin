#
# Anaquin - Sequin statistical analysis
#
# This R script was generated at %1%.
#
#    %2%
#

library(plyr)
library(Anaquin)

d1 <- read.csv('%3%/%4%', sep='\t')
d2 <- read.csv('%3%/%5%', sep='\t')

# Merge detected and known variants
data <- rbind.fill(d1, d2[d2$Label == 'FN',])

# Merge insertions and deletions
data$Mutation <- revalue(data$Mutation, c(Insertion = 'Indel', Deletion = 'Indel'))

# Construct unique identifiers
data$Unique <- paste(paste(data$Name, data$Pos, sep='_'), data$Type, sep='_')

data$Label <- revalue(data$Label, c("TP"="1", "FN"="1", 'FP'='0'))

plotROC(data$Unique, data$%6%, data$Mutation, data$Label, title='ROC (ranked by %6%)', refGroup=%8%, legTitle=%7%, 
        ords = c('0'='0', '1'='1'), xlab='FPR', ylab='TPR')
