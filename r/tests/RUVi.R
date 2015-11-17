#
#  Copyright (C) 2015 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

library(RUnit)
library(RUVSeq)

#
# Unit tests for RUVi normalization with TransQuin isoforms
#

.genes <- function()
{
    d <- c('/Users/tedwong/Dropbox/Sequins/Manuscripts/RNA/RUV/R1000/A1/genes.fpkm_tracking',
           '/Users/tedwong/Dropbox/Sequins/Manuscripts/RNA/RUV/R1000/A2/genes.fpkm_tracking',
           '/Users/tedwong/Dropbox/Sequins/Manuscripts/RNA/RUV/R1000/A3/genes.fpkm_tracking',
           '/Users/tedwong/Dropbox/Sequins/Manuscripts/RNA/RUV/R1000/B1/genes.fpkm_tracking',
           '/Users/tedwong/Dropbox/Sequins/Manuscripts/RNA/RUV/R1000/B2/genes.fpkm_tracking',
           '/Users/tedwong/Dropbox/Sequins/Manuscripts/RNA/RUV/R1000/B3/genes.fpkm_tracking')
    d
}

.isoforms <- function()
{
    d <- c('/Users/tedwong/Dropbox/Sequins/Manuscripts/RNA/RUV/R1000/A1/isoforms.fpkm_tracking',
           '/Users/tedwong/Dropbox/Sequins/Manuscripts/RNA/RUV/R1000/A2/isoforms.fpkm_tracking',
           '/Users/tedwong/Dropbox/Sequins/Manuscripts/RNA/RUV/R1000/A3/isoforms.fpkm_tracking',
           '/Users/tedwong/Dropbox/Sequins/Manuscripts/RNA/RUV/R1000/B1/isoforms.fpkm_tracking',
           '/Users/tedwong/Dropbox/Sequins/Manuscripts/RNA/RUV/R1000/B2/isoforms.fpkm_tracking',
           '/Users/tedwong/Dropbox/Sequins/Manuscripts/RNA/RUV/R1000/B3/isoforms.fpkm_tracking')
    d
}

testIsoforms <- function()
{
    i <- .isoforms()
    
    # Construct a count table at the isoform level
    d <- round(cufflink.isoforms(i[[1]], i[[2]], i[[3]], i[[4]], i[[5]], i[[6]]))

    # What does the RLE plot look like before normalization? In a real experiment, we would probably need a filter.
    plotRLE(d)

    # What does the PCA plot look like before normalization?
    plotPCA(d)
    
    # Let's normalize the counts by isoform FPKMs
    r <- TransNorm(d, level='isoforms')

    # What does the RLE plot look like after normalization?
    plotRLE(r$normalizedCounts)

    # What does the PCA plot look like after normalization?
    plotPCA(r$normalizedCounts)    
}

testGenes <- function()
{
    i <- .genes()

    # Construct a count table at the gene level
    d <- round(cufflink.genes(i[[1]], i[[2]], i[[3]], i[[4]], i[[5]], i[[6]]))
    
    # What does the RLE plot look like before normalization? In a real experiment, we would probably need a filter.
    plotRLE(d)
    
    # What does the PCA plot look like before normalization?
    plotPCA(d)
    
    # Let's normalize the counts by gene FPKMs
    r <- TransNorm(d, level='genes', method='neg')
    
    # What does the RLE plot look like after normalization?
    plotRLE(r$normalizedCounts)
    
    # What does the PCA plot look like after normalization?
    plotPCA(r$normalizedCounts)    
}
