#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

# Load reference sequins
data <- read.csv('%3%/%4%', row.names=1, sep='\t')

# Classify the sequins against the negative controls
data$label <- classifyByRefRatio(inputs=abs(data$ExpLFC), refRatio=0)

# Create Anaquin data set
data <- createAnaquinData(names=row.names(data), input=data$ExpLFC, measured=data$ObsLFC, score=1-data$Pval, qval=data$Pval, label=data$label)

# Title of the ROC plot
title <- 'ROC Plot'

plotROC(data, title=title, refRats=0)
