#
#  Copyright (C) 2015 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
#

library(RUnit)
library(RUVSeq)

testRMXA <- function()
{
    setwd('/Users/tedwong/Sources/QA/r/tests')    
    d  <- read.csv('data/data.csv', row.names=1)
    colnames(d) <- c('A1', 'A2', 'A3', 'B1', 'B2', 'B3')

    filter <- apply(d, 1, function(x) length(x[x>5])>=2)
    d <- d[filter,]

    m <- loadMixture()
    detected <- rownames(d) %in% m$genes$ID
    spikes   <- rownames(d[detected,])

    r1 <- RUVg(as.matrix(d), spikes, k=1)
    r2 <- TransNorm(d)
    
    checkTrue(identical(r1,r2))
}

testRMXA_1 <- function()
{
    setwd('/Users/tedwong/Sources/QA/r/tests')    
    d  <- read.csv('data/data.csv', row.names=1)
    colnames(d) <- c('A1', 'A2', 'A3', 'B1', 'B2', 'B3')
    
    filter <- apply(d, 1, function(x) length(x[x>5])>=2)
    d <- d[filter,]
    
    x <- c('R2_54', 'R1_23', 'R2_7', 'R1_71', 'R2_152', 'R2_18', 'R1_81', 'R2_117', 'R2_45')
    
    m <- loadMixture()
    #m <- m$genes[m$genes$Fold==1,]
    m <- m$genes[m$genes$ID %in% x,]
    
    detected <- rownames(d) %in% m$ID
    spikes   <- rownames(d[detected,])
    
    r1 <- RUVg(as.matrix(d), spikes, k=2)
    r2 <- TransNorm(d)
    
    checkTrue(identical(r1,r2))
}
