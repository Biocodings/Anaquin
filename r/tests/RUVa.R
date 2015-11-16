#
#  Copyright (C) 2015 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
#

library(RUnit)
library(RUVSeq)

#
# Unit tests for RUVa normalization with TransQuin sequins.
#

.data <- function()
{
    d  <- read.csv('/Users/tedwong/Sources/QA/r/tests/data/data.csv', row.names=1)
    d
}

#
# Test with all sequins, note that this is generally not a good idea but we're only concern the result here
#
testRMXA_1 <- function()
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
    r2 <- TransNorm(d, method='all')

    checkTrue(identical(r1,r2))
}

#
# Repeat the previous test but with a subset of sequins
#
testRMXA_2 <- function()
{
    setwd('/Users/tedwong/Sources/QA/r/tests')    
    d  <- read.csv('data/data.csv', row.names=1)
    colnames(d) <- c('A1', 'A2', 'A3', 'B1', 'B2', 'B3')
    
    filter <- apply(d, 1, function(x) length(x[x>5])>=2)
    d <- d[filter,]
    
    x <- c('R2_54', 'R1_23', 'R2_7', 'R1_71', 'R2_152', 'R2_18', 'R1_81', 'R2_117', 'R2_45')
    
    m <- loadMixture()
    m <- m$genes[m$genes$ID %in% x,]
    
    detected <- rownames(d) %in% m$ID
    spikes   <- rownames(d[detected,])
    
    r1 <- RUVg(as.matrix(d), spikes, k=1)
    r2 <- TransNorm(d, spikes=spikes)
    
    checkTrue(identical(r1,r2))
}

#
# What happens to a perfect experiement where all the negative control sequins really are truly unaffected?
#
testRMXA_3 <- function()
{
    d <- .data()
    m <- loadMixture()

    detected <- rownames(d) %in% m$genes$ID

    # This is used to generate pseduo counts    
    i <- 1
    
    for (id in row.names(d[detected,]))
    {
        d[id,]$A1 <-  (10 * i) + rnorm(1, 5, 2)
        d[id,]$A2 <-  (10 * i) + rnorm(1, 5, 2)
        d[id,]$A3 <-  (10 * i) + rnorm(1, 5, 2)
        d[id,]$B1 <-  (10 * i) + rnorm(1, 5, 2)
        d[id,]$B2 <-  (10 * i) + rnorm(1, 5, 2)
        d[id,]$B3 <-  (10 * i) + rnorm(1, 5, 2)
        i <- i + 1
    }
    
    r <- TransNorm(d, method='all')
}

#
# What happens to a perfect experiement where all the negative control sequins have been biased two-times?
#
testRMXA_4 <- function()
{
    d <- .data()
    m <- loadMixture()
    
    # This is used to generate pseduo counts    
    i <- 1
    
    for (id in row.names(d[rownames(d) %in% m$genes$ID,]))
    {
        d[id,]$A1 <- (10 * i) + rnorm(1, 5, 2)
        d[id,]$A2 <- (10 * i) + rnorm(1, 5, 2)
        d[id,]$A3 <- (10 * i) + rnorm(1, 5, 2)
        d[id,]$B1 <- 2 * d[id,]$A1
        d[id,]$B2 <- 2 * d[id,]$A2
        d[id,]$B3 <- 2 * d[id,]$A3
        i <- i + 1
    }

    r <- TransNorm(d, method='all')
}
