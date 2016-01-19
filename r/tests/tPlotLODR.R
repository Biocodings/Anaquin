#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

library('RUnit')
library('Anaquin')

plotForGenes <- function()
{
    d <- read.csv('/Users/tedwong/Desktop/LODR_genes_pval_TED.csv', row.names=1)
    d <- d[!is.na(d$pvalue),]
    
    data <- transQuin(seqs = row.names(d), counts = d$baseMean, pval = d$pvalue, ratio = d$expected.LFC)

    plotLODR(data, shouldTable=FALSE, lvl='gene', shouldBand=TRUE)
}

plotForExons <- function()
{
    #
    # Create an LODR plot for exons. The data set is the same data set used in the sequins paper.
    #
    d <- read.csv('/Users/tedwong/Desktop/LODR_exons_TED.csv', row.names=1)
    d$LFC <- round(d$LFC)
    
    d <- d[!is.na(d$pvalue),]

    #
    # The following LFC can't be fitted by local regression
    #
    #     - (-8,8): 3 counts
    #     - (-6,6): 2 counts
    #
    
    d <- d[d$LFC!=8 & d$LFC!=-8,]
    d <- d[d$LFC!=6 & d$LFC!=-6,]
    d <- d[d$LFC!=7 & d$LFC!=-7,]
    d <- d[d$LFC!=9 & d$LFC!=-9,]    

    #
    # It's hard to construct a LODR plot with zero probabilities...
    #

    #m <- mean(d[(d$LFC==-5 | d$LFC==5) & (d$pvalue!=0),]$pvalue)
    #d[(d$LFC==-5 | d$LFC==5) & (d$pvalue==0),]$pvalue <- 0.001
    
    d[d$pvalue==0,]$pvalue <- 1e-100
    

#    d <- d[d$pvalue>0,]
    d <- d[d$pvalue>1e-200,]

    data <- anaquin(seqs=row.names(d), counts=d$exonBaseMean, pval=d$pvalue, ratio=d$LFC)
    plotLODR(data, shouldTable=FALSE)
}

plotForExons()
plotForGenes()
