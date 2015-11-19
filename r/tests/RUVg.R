#
#  Copyright (C) 2015 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
#

library(RUnit)
library(RUVSeq)

#
# Unit tests for RUVg normalization with TransQuin genes
#

.data <- function()
{
    d  <- read.csv('/Users/tedwong/Dropbox/Sequins/Manuscripts/RNA/RUV/Real/genes.csv', row.names=1)
    colnames(d) <- c('A1', 'A2', 'A3', 'B1', 'B2', 'B3')
    d
}

.sequins <- function()
{
    d <- .filter(.data(), loadMixture())
    d
}

#
# Demonstrate how RUVg can be used to normalize samples with sequins, with the real data
#
caseStudy <- function()
{
    # Load the data, including real gene and sequins
    d <- .data()

    # Filter by 5 reads
    d <- d[apply(d, 1, function(x) length(x[x>5])>=2),]

    # Negative control genes
    spikes <- c('R1_21', 'R1_23', 'R1_71', 'R1_81', 'R2_117', 'R2_140', 'R2_152', 'R2_18','R2_20','R2_45','R2_54','R2_65','R2_7','R2_71')
    #spikes <- c('R1_21')
    
    #
    # Let's try to reuse the original package. To make things simpler, we'll prefer SeqExpressionSet.
    #

    x2 <- as.factor(rep(c('A', 'B'), each=3))
    s2 <- newSeqExpressionSet(as.matrix(d), phenoData=data.frame(x2, data.frame(x2, row.names=colnames(d))))

    library(RColorBrewer)
    colors <- brewer.pal(3, 'Set2')
    EDASeq::plotRLE(s2, outline=FALSE, ylim=c(-4, 4), col=colors[x])
    EDASeq::plotPCA(s2, col=colors[x], cex=1.2)

    r2 <- RUVSeq::RUVg(s2, rownames(d[rownames(d) %in% spikes,]), k=1)
    EDASeq::plotRLE(r2, outline=FALSE, ylim=c(-4, 4), col=colors[x])
    EDASeq::plotPCA(r2, col=colors[x], cex=1.2)

    #
    # Let's repeat the analysis with our RUV implementation
    #
        
    r1 <- TransNorm(d, level='genes', method='neg', spikes=spikes)
    x1  <- as.factor(rep(c('A', 'B'), each=3))
    s1 <- newSeqExpressionSet(as.matrix(r1$normalizedCounts), phenoData=data.frame(x1, data.frame(x1, row.names=colnames(d))))

    EDASeq::plotRLE(s1, outline=FALSE, ylim=c(-4, 4), col=colors[x])
    EDASeq::plotPCA(s1, col=colors[x], cex=1.2)
}

#
# Test with the negative controls at the gene level. No filter is needed since only sequins are involved.
#
test_1 <- function()
{
    # Load Simon's data, obviously it's unnormalized
    d <- .sequins()

    # What does the RLE plot look like before normalization? In a real experiment, we would probably need a filter.
    plotRLE(d)
    
    # What does the PCA plot look like before normalization?
    plotPCA(d)
    
    # Let's normalize the counts by genes
    r <- TransNorm(d, level='genes')
    
    # What does the RLE plot look like after normalization?
    plotRLE(r$normalizedCounts)
    
    # What does the PCA plot look like after normalization?
    plotPCA(r$normalizedCounts)
}

#
# Test with all sequins, note that this is generally not a good idea but we're only concern the result here
#
test_2 <- function()
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
test_3 <- function()
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
test_4 <- function()
{
    set.seed(1234)

    d <- .data()
    m <- loadMixture()

    # This is used to generate pseduo counts    
    i <- 1
    
    for (id in row.names(d[rownames(d) %in% m$genes$ID,]))
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

    checkTrue(sum(c(r$normalizedCounts['R2_76',])) == 4595)
}

#
# What happens to a perfect experiement where all the negative control sequins have been biased two-times?
#
test_5 <- function()
{
    set.seed(1234)

    d <- .data()
    m <- loadMixture()
    
    # This is used to generate pseduo counts    
    i <- 1
    
    #
    # Let's pretend we have a library bias (eg: library size) that makes the second mixture have twice
    # number of counts. This is an unwanted variation, our RUV method should be able to normalize the
    # counts back to what it should have been if there was no bias. For example, the fold-change of two
    # could make differentially expressed, when it shouldn't.
    #

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
    checkTrue(all(c(r$normalizedCounts['R2_76',]) == c(1086, 1088, 1081, 1086, 1088, 1081)))
}
