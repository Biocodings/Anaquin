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
    # Create an MA plot for exons. The data set is the one used in the sequins paper. Contributed
    # by Simon Hardwick at Garvan Institute.
    #
    
    #
    # 1. Read the counts for sequin exon bins
    #
    
    d <- read.csv('/Users/tedwong/Desktop/LODR_exons_TED.csv', row.names=1)
    colnames(d) <- c('Mean', 'PValue', 'QValue', 'B1', 'B2', 'B3', 'A1', 'A2', 'A3', 'logFold')
    d$logFold <- round(d$logFold)
    seqs <- data.frame(A1=seqs$A1, A2=seqs$A2, A3=seqs$A3, B1=seqs$B1, B2=seqs$B2, B3=seqs$B3, ratio=seqs$logFold)
    row.names(seqs) <- row.names(d)

    #
    # 2. Read the counts for K_562 exon bins
    #
    
    K_562 <- read.csv('/Users/tedwong/Desktop/Exon_K562_Counts.csv', row.names=1)
    K_562$ratio <- NA

    stopifnot(sum(is.na(seqs$ratio))   == 0)
    stopifnot(sum(!is.na(K_562$ratio)) == 0)    
        
    data <- rbind(seqs, K_562)
    
    r <- plotMA(data=data, metrs='exon', shouldEndo=TRUE, shouldError=TRUE)

    checkEquals(r$xname, 'Log2 Average of Normalized Counts')
    checkEquals(r$yname, 'Log2 Ratio of Normalized Counts')
}

plotForGenes <- function()
{
    d <- read.csv('tests/data/counts.txt', row.names=1)
    row.names(d) <- d$Feature
    d <- d[,-1]

    r <- plotMA(d, metrs='gene', shouldError=TRUE, shouldEndo=TRUE)
    
    checkEquals(r$xname, 'Log2 Average of Normalized Counts')
    checkEquals(r$yname, 'Log2 Ratio of Normalized Counts')
}

plotForExons()
plotForGenes()
