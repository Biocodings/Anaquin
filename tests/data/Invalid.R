#
# Anaquin - Sequin statistical analysis. Version 1.1.1.
#
# This R script was generated at 01-08-2016 12:55:59.
#
#    anaquin RnaExpression -o ABCD -m MRN027_v001.csv -method gene -ufiles A1/G/transcripts.gtf 
#

library('Anaquin')

# Load reference sequins
data <- read.csv('/data/Pipeline/RnaQuin_2/ABCD/RnaExpression_sequins.csv', row.names=1, sep='\t')

title <- 'Gene Expression'
xlab  <- 'Input Concentration (log2)'
ylab  <- 'FPKM (log2)'

# Expected input concentration (x-axis)
input <- log2(data$InputConcent)

# Measured expression (y-axis)
measured <- log2(data$Observed)

# Create Anaquin data for plotScatter
anaquin <- createAnaquinData(names=row.names(data), input=input, measured=measured)

plotScatter(anaquin, title=title, xlab=xlab, ylab=ylab, showLOQ=TRUE)$
