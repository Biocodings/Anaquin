#
#  Copyright (C) 2015 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
#

library(RUnit)
library(RUVSeq)

#
# Unit tests for RUVe normalization with TransQuin exons
#

.data <- function()
{
    d  <- read.csv('/Users/tedwong/Sources/QA/r/tests/data/exons.csv', row.names=1)
    d
}

testExons_1 <- function()
{
    d <- .data()   
    
    
    
    
    colnames(d) <- c('A1', 'A2', 'A3', 'B1', 'B2', 'B3')

    filter <- apply(d, 1, function(x) length(x[x>5])>=2)
    d <- d[filter,]

    m <- loadMixture()
    detected <- rownames(d) %in% m$genes$ID
    spikes   <- rownames(d[detected,])

    r1 <- RUVg(as.matrix(d), spikes, k=1)
    r2 <- TransNorm(d, method='all')

    checkTrue(identical(r1,r2))
}