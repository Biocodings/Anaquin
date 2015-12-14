#
#  Copyright (C) 2015 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

library(DESeq2)

# Load the count matrix
d <- read.csv('data.csv', row.names=1)

meta <- read.csv(file.path('tutorial.csv'), row.names=1)

se <- DESeqDataSetFromMatrix(d, DataFrame(meta), design=~Sample)

dds <- DESeqDataSet(se, design = ~Sample)
dds <- DESeq(dds)
res <- results(dds, contrast=c("Sample","G_RMXB","K_RMXA"))

#
# Now, try to gain some insights with Anaquin. 
#

r <- .analyzeDESeq2(res, loadMixture())

plotROC(r$data)





