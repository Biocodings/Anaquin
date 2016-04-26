#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

library('RUnit')
library('DESeq2')
library('Anaquin')

testTR04449 <- function()
{
    meta <- read.csv(file.path('tests/data/meta.csv'), row.names=1)

    data <- read.csv('tests/data/counts.csv', row.names=1) 
    se <- DESeqDataSetFromMatrix(data, DataFrame(meta), design=~Sample)
    
    dds <- DESeqDataSet(se, design = ~Sample)
    dds <- DESeq(dds)
    r   <- results(dds, contrast=c("Sample", "G_RMXB", "K_RMXA"))
    
    data <- TransQuin(mix='AB')
    TransDiff(data, r)
}

testTR04449()