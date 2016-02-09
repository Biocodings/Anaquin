#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
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
