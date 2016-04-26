#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

library('RUnit')
library('edgeR')
library('Anaquin')
library('GenomicFeatures')

testTR04449_edgeR <- function()
{
    data  <- read.csv('tests/data/counts.csv', row.names=1)
    group <- factor(c(1,1,1,2,2,2))
    
    y <- DGEList(counts=assay(se), group=group)
    y <- calcNormFactors(y)
    y <- estimateCommonDisp(y)
    y <- estimateTagwiseDisp(y)
    r <- exactTest(y)
    
    data <- TransQuin(mix='AB')
    TransDiff(data, r)
}

testTR04449_edgeR()