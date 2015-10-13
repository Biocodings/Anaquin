#
#  Copyright (C) 2015 - Garvan Institute of Medical Research
#
#  Written by Ted Wong, Bioinformatic Engineer at Garvan Institute.
#

#
# Tools for visualizing synthetic chromosome 
#
#    - Density plot for coverage
#

library('Sushi')
library(rtracklayer)
library(GenomicRanges)

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
        plotBedgraph(src, 'chrT', start(seq), end(seq))
    }
}

plotDensity('/Users/tedwong/Sources/QA/output/VarCoverage_summary.bedgraph', '/Users/tedwong/Sources/QA/data/var/AVA017.v032.bed')




