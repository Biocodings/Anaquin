#
#  Copyright (C) 2015 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
#

library(RUnit)
library(RUVSeq)

testRMXAv2 <- function()
{
    #
    # We'll manually manipulate the data and then compare the result with Anaquin.
    #

    d  <- read.csv('data/data.csv', row.names=1)
    colnames(d) <- c('A1', 'A2', 'A3', 'B1', 'B2', 'B3')

    filter <- apply(d, 1, function(x) length(x[x>5])>=2)
    d <- d[filter,]

    m  <- loadMixture()
    detected <- rownames(d) %in% m$genes$ID
    spikes   <- rownames(d[detected,])

    r1 <- RUVg(as.matrix(d), spikes, k=1)
    r2 <- TransNorm(d)
    
    checkTrue(identical(r1,r2))
}


