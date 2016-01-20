#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

library('RUnit')
library('Anaquin')

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

plotForGenes <- function()
{
    # Read data file for gencode
    gens <- read.csv('tests/data/K_562/DESeq2_gencode_results.csv', row.names=1)
    
    # Read data file for sequins
    seqs <- read.csv('tests/data/K_562/LODR_genes_TED_20.01.16.csv', row.names=1)
    
    # There is no expected LFC for gencode but we'll need it to do a column bind
    gens$expected.LFC <- NA
    
    data <- rbind(seqs[,c(1,2,3,5,6),], gens[,c(1,2,3,5,7)])
    stopifnot((nrow(seqs) + nrow(gens)) == nrow(data))
    
    # Create a TransQuin data set for Anaquin
    data <- TransQuin(seqs=row.names(data), baseMean=data$baseMean, log2FoldChange=data$log2FoldChange, pvalue=data$pvalue, expected.LFC=data$expected.LFC)

    plotLODR(data, choseFDR=0.1, shouldTable=FALSE, lvl='gene', shouldBand=FALSE, yBreaks=c(1e-300, 1e-200, 1e-100, 1e-10, 1.00), locBand='local')
}

#plotForExons()
plotForGenes()
