#
#  Copyright (C) 2015 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

library(RUnit)
library(RUVSeq)
library(DEXSeq)

#
# Unit tests for RUVe normalization with TransQuin exons
#

.files <- function()
{
    d <- c('~/Dropbox/Sequins/Manuscripts/RNA/RUV/Real/A1.txt',
           '~/Dropbox/Sequins/Manuscripts/RNA/RUV/Real/A2.txt',
           '~/Dropbox/Sequins/Manuscripts/RNA/RUV/Real/A3.txt',
           '~/Dropbox/Sequins/Manuscripts/RNA/RUV/Real/B1.txt',
           '~/Dropbox/Sequins/Manuscripts/RNA/RUV/Real/B2.txt',
           '~/Dropbox/Sequins/Manuscripts/RNA/RUV/Real/B3.txt')
    d
}

.sampleData <- function()
{
    d <- data.frame(row.names = c("A1", "A2", "A3", "B1", "B2", "B3"), condition = c("A", "A", "A", "B", "B", "B"))
    d
}

.flattenGTF <- function()
{
    r <- '~/Dropbox/Sequins/Manuscripts/RNA/RUV/Real/combined.flattened.gtf'
    r
}

#
# This test demostrates a simple usage of normalization by exons
#

testExons <- function()
{
    d <- DEXSeqDataSetFromHTSeq(.files(), sampleData=.sampleData(), design=~sample+exon+condition:exon, flattenedfile=.flattenGTF())
    f <- counts(d)
    f <- f[,c(1:6)]
    f <- f[apply(f, 1, function(x) length(x[x>5])>=2),]
    colnames(f) <- c('A1', 'A2', 'A3', 'B1', 'B2', 'B3')

    x <- as.factor(rep(c('A', 'B'), each=3))
    s <- newSeqExpressionSet(as.matrix(f), phenoData=data.frame(x2, data.frame(x, row.names=colnames(f))))
    
    library(RColorBrewer)
    colors <- brewer.pal(3, 'Set2')
    EDASeq::plotRLE(s, outline=FALSE, ylim=c(-4, 4), col=colors[x])
    EDASeq::plotPCA(s, col=colors[x], cex=1.2)

    # Negative control exon bins
    spikes <- negativeExonBins(loadMixture())

    r2 <- RUVSeq::RUVg(s, rownames(f[rownames(f) %in% row.names(spikes),]), k=1)

    EDASeq::plotRLE(r2, outline=FALSE, ylim=c(-4, 4), col=colors[x])
    EDASeq::plotPCA(r2, col=colors[x], cex=1.2)
    
    # Let's normalize the counts by exon bins
    #r <- TransNorm(f, level='exons')

    #EDASeq::plotRLE(r$normalizedCounts)
    #EDASeq::plotPCA(r$normalizedCounts)
}
