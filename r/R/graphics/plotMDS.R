#
#  Copyright (C) 2015 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

#
# ----------------------- ROC Plot -----------------------
#

plotROC <- function(r, labels)
{
    d <- r$data
    d <- d[!is.nan(d$measured),]
    d <- d[d$class == 'TP' | d$class == 'FP',]

    # 
    # From the package's reference manual:
    #
    # ... labels should be supplied as ordered factor(s), the lower level corresponding to the negative class, the upper level
    #     to the positive class ...
    #
    d$label <- ifelse(d$class == 'TP', 2, 1)
    
    d$scores <- 1 - d$pval

    require(ROCR)

    scores <- c(0.999985006663244,0.11848973121444,0.998317078130317,0.209432682723151,0.991844504940773,0.296175974698939,0.890247946913349,0.67291092334451,0.324408809166977,0.103365779196362,0.252794403922691,0.600822276080659,0.045599350395611,0.99997718730894,0.543686960607272,0.806223370592858,0.999698543215618,0.995835949956063,0.183529103779826,0.99986388387351,0.130236027425796,0.999552354071581,0.999917043997976,0.999911464673142,0.99985115274586,0.976245886879838,0.341646255148427,0.825269252328571,0.14267747036437,0.987840216112876,0.033142220590023)
    labels <- c(4,1,4,1,4,1,4,1,1,1,1,1,1,4,1,4,4,4,1,4,1,4,4,4,4,4,1,4,1,4,1)
    
    pred <- prediction(scores, labels, label.ordering=c(1,4))
    #pred <- prediction(d$scores, d$label, label.ordering=c(1,2))
    perf <- performance(pred, "tpr","fpr")
    
#    plot(perf, colorize=TRUE, cex.lab=2, cex.main=2, lwd=10)    
    plot(perf)
}

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