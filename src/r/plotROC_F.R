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
data <- FusQuin(seqs=data$Sequin, label=data$Label, pos1=data$Position_1, pos2=data$Position_2)

#plotROC.FusQuin(data)