#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

library('RUnit')
library('Anaquin')

DESeq2 <- function()
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

    xBreaks <- c(1, 10, 100, 1000, 10000)
    xLabels <- c('1e+00', '1e+01', '1e+02', '1e+03', '1e+04')
    yBreaks <- c(1e-300, 1e-200, 1e-100, 1e-10, 1.00)

    r <- plotLODR(data, shouldBand=FALSE, locBand='pred', choseFDR=0.1, shouldTable=FALSE, lvl='gene', yBreaks=yBreaks, xBreaks=xBreaks, xLabels=xLabels)
}

DESeq2()