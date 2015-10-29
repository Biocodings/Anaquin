#
#  Copyright (C) 2015 - Garvan Institute of Medical Research
#
#  Written by Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
#

#
# Tools for visualizing synthetic chromosome 
#
#    - Scatter plot for concentration vs coverage
#    - Density plot for coverage
#

library('Sushi')
library('ggplot2')
library(rtracklayer)
library(GenomicRanges)
library(plyr)

#
# Scatter plot is the most common data visualization tool in Anaquin. It plots the expected concentration
# defined by a mixture file with the measured coverage.
#

plotScatter <- function(x, y, ids, isLog=FALSE)
{
    if (!isLog)
    {
        d <- data.frame(x=log2(x), y=log2(y), ids=ids)
    }
    else
    {
        d <- data.frame(x=x, y=y, ids=ids)
    }

    p <- ggplot(data = d, aes(x = x, y = y))
    p <- p + xlab('Expected log2 fold change of mixture A and B')
    p <- p + ylab('Measured log2 fold change of mixture A and B')
    p <- p + geom_point()
    p <- p + ggtitle('')
    p <- p + xlim(min(d$x)-2, max(d$x)+2)
    p <- p + ylim(min(d$y)-2, max(d$y)+2)
    p <- p + geom_smooth(method = 'lm', formula = y ~ x)
    #p + geom_text(x = 0, y = max(d$y), label = lm_eqn(d), parse = TRUE)
    print(p)
}

plotDensity <- function(src, ref)
{
    # Read the source file
    src <- read.csv(src, header=FALSE, sep='\t')
    
    colnames(src) <- c('chrom', 'start',  'end',  'value')
    
    # Read the reference annotation
    ref <- import.bed(con=ref, asRangedData=F)
    
    # Loop over all sequins defined in the reference...
    for (i in 1:length(ref))
    {
        seq <- ref[i,]

        # Construct a density plot for the sequin
        plotBedgraph(src,
                     'chrT',
                     start(seq),
                     end(seq),
                     transparency=0.50,
                     color='#ADC2E6',
                     xlab=seq$name,
                     ylab='Coverage')

        ticks  <- 5
        range  <- c(min(src[src$start >= start(seq) & src$end <= end(seq),]$value),
                    max(src[src$start >= start(seq) & src$end <= end(seq),]$value))
        scaled <- range / ticks
        scaled <- round_any(scaled, 100, f = ceiling)

        axis(side=2, at=seq(0, ticks * scaled[2], by=scaled[2]))
    }
}

#plotDensity('/Users/tedwong/Sources/QA/output/VarCoverage_summary.bedgraph', '/Users/tedwong/Sources/QA/data/var/AVA017.v032.bed')
