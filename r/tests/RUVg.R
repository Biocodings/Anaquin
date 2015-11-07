#
#  Copyright (C) 2015 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
#

#
# This file defines tests for the RUVg normalization.
#

library(RUVSeq)
library("RUnit")

library(RColorBrewer)
colors <- brewer.pal(3, "Set2")

testRMXAv2 <- function()
{
    d <- read.csv('/Users/tedwong/Desktop/data.csv', row.names=1)
    colnames(d) <- c('A1', 'A2', 'A3', 'B1', 'B2', 'B3')

    filter <- apply(d, 1, function(x) length(x[x>5])>=2)
    d <- d[filter,]

    x <- as.factor(rep(c("MixA", "MixB"), each=3))
    before <- newSeqExpressionSet(as.matrix(d), phenoData=data.frame(x, row.names=colnames(d)))
    plotPCA(before)
    plotRLE(before, outline=FALSE, col=colors[x])
    
    r <- TransNorm(d)
    
    after <- newSeqExpressionSet(as.matrix(r$normalizedCounts), phenoData=data.frame(x, row.names=colnames(r$normalizedCounts)))
    plotPCA(after)
    plotRLE(after, outline=FALSE, col=colors[x])
}
