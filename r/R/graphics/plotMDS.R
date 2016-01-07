#
#  Copyright (C) 2015 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

#
# ----------------------- PCA Plot -----------------------
#
# Create a scatter plot for the first two principal components for RUVg.
#

plotPCA <- function(m)
{
    x <- as.factor(rep(c("MixA", "MixB"), each=3))
    s <- newSeqExpressionSet(as.matrix(m), phenoData=data.frame(x, row.names=colnames(m)))
    EDASeq::plotPCA(s, outline=FALSE, col=colors[x])
}

#
# ----------------------- RLE Plot -----------------------
#
# Create a plot for relative likehihood expression for RUVg.
#

plotRLE <- function(m)
{
    x <- as.factor(rep(c("MixA", "MixB"), each=3))
    s <- newSeqExpressionSet(as.matrix(m), phenoData=data.frame(x, row.names=colnames(m)))    
    EDASeq::plotRLE(s, outline=FALSE, col=colors[x])
}

#
# ----------------------- Density Plot -----------------------
#
# Density plot tabulates the coverage across sequins.
#plotDensity('/Users/tedwong/Sources/QA/TransCoverage_chrT.bedgraph', '/Users/tedwong/Sources/QA/data/trans/ATR002.v032.bed', seqIDs=c('R2_14_1'))
#plotDensity('/Users/tedwong/Sources/QA/TransCoverage_chrT.bedgraph', '/Users/tedwong/Sources/QA/data/trans/ATR002.v032.bed')

plotDensity <- function(src, ref, seqIDs = NULL, minBase = NULL, maxBase = NULL)
{
    require(Sushi)
    
    # Read the source file
    src <- read.csv(src, header=FALSE, sep='\t')
    
    colnames(src) <- c('chrom', 'start',  'end',  'value')
    
    # Read the reference annotation
    ref <- import.bed(con=ref)
    
    # Loop over all sequins defined...
    for (i in 1:length(ref))
    {
        seq_ <- ref[i,]
        
        if (is.null(seqIDs) || seq_$name %in% seqIDs)
        {
            startL <- start(seq_)
            endL   <- end(seq_)
            
            if ((is.null(minBase) || startL >= minBase) && (is.null(maxBase) || endL <= maxBase))
            {
                Sushi::plotBedgraph(src, 'chrT', start(seq_), end(seq_), transparency=0.50, color='#ADC2E6', xlab=seq_$name, ylab='Coverage')
                #Sushi::plotBedgraph(src, 'chrT', 6535259, 6536017, transparency=0.50, color='#ADC2E6', xlab=seq_$name, ylab='Coverage')
                
                ticks  <- 5
                range  <- c(min(src[src$start >= start(seq_) & src$end <= end(seq_),]$value),
                            max(src[src$start >= start(seq_) & src$end <= end(seq_),]$value))
                scaled <- range / ticks
                scaled <- round_any(scaled, 100, f = ceiling)
                
                axis(side=2, at=seq(0, ticks * scaled[2], by=scaled[2]))
            }
        }
    }
}