#
# Anaquin - Sequin statistical analysis
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

data <- read.csv('%3%/%4%', row.names=1, sep='\t')

# Classify sequins against the negative controls (LFC 0)
data$label <- ifelse(abs(data$ExpLFC) <= 0, '0', '1')

title <- 'ROC Plot'

# Create Anaquin data for PlotROC
anaquin <- AnaquinData(analysis='PlotROC', seqs=row.names(data), ratio=data$ExpLFC, measured=data$ObsLFC, score=1-data$Pval, label=data$label)

plotROC(anaquin, title=title, refRats=0)
