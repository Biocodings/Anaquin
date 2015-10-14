#
#  Copyright (C) 2015 - Garvan Institute of Medical Research
#
#  Written by Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
#

#
# Tools for visualizing synthetic chromosome 
#
#    - Density plot for coverage
#

library('Sushi')
library(rtracklayer)
library(GenomicRanges)
library(plyr)

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

plotDensity('/Users/tedwong/Sources/QA/output/VarCoverage_summary.bedgraph', '/Users/tedwong/Sources/QA/data/var/AVA017.v032.bed')
