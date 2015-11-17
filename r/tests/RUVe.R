#
#  Copyright (C) 2015 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
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

.files <- function()
{
    d <- c('/Users/tedwong/Sources/QA/r/tests/data/simon_experiment/K_RMXA1v2.DEXSeq.counts.txt',
           '/Users/tedwong/Sources/QA/r/tests/data/simon_experiment/K_RMXA2v2.DEXSeq.counts.txt',
           '/Users/tedwong/Sources/QA/r/tests/data/simon_experiment/K_RMXA3v2.DEXSeq.counts.txt',
           '/Users/tedwong/Sources/QA/r/tests/data/simon_experiment/G_RMXB1v2.DEXSeq.counts.txt',
           '/Users/tedwong/Sources/QA/r/tests/data/simon_experiment/G_RMXB2v2.DEXSeq.counts.txt',
           '/Users/tedwong/Sources/QA/r/tests/data/simon_experiment/G_RMXB3v2.DEXSeq.counts.txt')
    d
}

.sampleData <- function()
{
    d <- data.frame(row.names = c("A1", "A2", "A3", "B1", "B2", "B3"), condition = c("A", "A", "A", "B", "B", "B"))
    d
}

.flattenGTF <- function()
{
    r <- '~/Sources/QA/data/trans/ATR001.v032.flatten.gtf'
    r
}

#
# This test demostrates a simple usage of normalization by exons
#

testExons_1 <- function()
{
    # Load Simon's data, obviously it's unnormalized
    d <- DEXSeqDataSetFromHTSeq(.files(), sampleData=.sampleData(), design=~sample+exon+condition:exon, flattenedfile=.flattenGTF())
    
    f <- counts(d)
    f <- f[,c(1:6)]
    colnames(f) <- c('A1', 'A2', 'A3', 'B1', 'B2', 'B3')

    # What does the RLE plot look like before normalization? In a real experiment, we would probably need a filter.
    plotRLE(f)
    
    # What does the PCA plot look like before normalization?
    plotPCA(f)

    # Let's normalize the counts by exon bins
    r <- TransNorm(f, level='exons')

    # What does the RLE plot look like after normalization?
    plotRLE(r$normalizedCounts)
    
    # What does the PCA plot look like after normalization?
    plotPCA(r$normalizedCounts)
}

