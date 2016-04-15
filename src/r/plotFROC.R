#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

data <- read.csv('/Users/tedwong/Sources/QA/output/FusionDiscover_labels.csv', sep='\t', stringsAsFactors=FALSE)

# Assign unique names to the false-positives
data[data$Label=='FP',]$Sequin <- c(1:nrow(data[data$Label=='FP',]))

# Create a FusQuin data set for Anaquin
data <- FusQuin(seqs=data$seq, label=data$label, pos1=data$pos1, pos2=data$pos2)

#plotROC.FusQuin(data)