#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

library('RUnit')
library('Anaquin')

plotForCustoms <- function()
{
    d <- read.csv('/Users/tedwong/Desktop/LODR_exons_TED.csv', row.names=1)
    colnames(d) <- c('Mean', 'PValue', 'QValue', 'B1', 'B2', 'B3', 'A1', 'A2', 'A3', 'logFold')
    d$logFold <- round(d$logFold)
        
    data <- data.frame(A1=d$A1, A2=d$A2, A3=d$A3, B1=d$B1, B2=d$B2, B3=d$B3, ratio=d$logFold)
    row.names(data) <- row.names(d)

    r <- plotMA(data, metrs='custom')

    checkEquals(r$xname, 'Log2 Average of Normalized Counts')
    checkEquals(r$yname, 'Log2 Ratio of Normalized Counts')
}

plotForGenes <- function()
{
    d <- read.csv('tests/data/counts.txt', row.names=1)
    row.names(d) <- d$Feature
    d <- d[,-1]
    
    r <- plotMA(d, metrs='gene', shouldEndo=TRUE)
    
    checkEquals(r$xname, 'Log2 Average of Normalized Counts')
    checkEquals(r$yname, 'Log2 Ratio of Normalized Counts')
}

plotForGenes()
plotForCustoms()
