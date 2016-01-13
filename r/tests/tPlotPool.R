#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

library('RUnit')
library('Anaquin')

testMinorMajor <- function()
{
    #
    # K_562 experiment for splicing minor/major at the gene level
    #

    data <- read.csv('/Users/tedwong/Desktop/K_RMXA_minIsoform_full_3reps_TED.csv', row.names=1)
    data <- data[-c(1),]
    colnames(data) <- c('X', 'A1', 'A2', 'A3')
    
    data$X  <- as.numeric(as.character(data$X))
    data$A1 <- as.numeric(as.character(data$A1))
    data$A2 <- as.numeric(as.character(data$A2))
    data$A3 <- as.numeric(as.character(data$A3))

    plotPool(data, cname='Ratio', xname='Expected log2 fold change of minor/major', yname='Measured log2 fold change of minor/major')
}

testMinorMajor()