#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

library('RUnit')
library('Anaquin')

testGenes <- function()
{
    data <- read.csv('tests/data/K_562/K_RMXA_minIsoform_full_3reps_TED.csv', row.names=1)
    data <- data[-c(1),]
    colnames(data) <- c('X', 'A1', 'A2', 'A3')
    
    data <- data[data$A1!='#DIV/0!',]
    data <- data[data$A2!='#DIV/0!',]
    data <- data[data$A3!='#DIV/0!',]
    
    data$X  <- as.numeric(as.character(data$X))
    data$A1 <- as.numeric(as.character(data$A1))
    data$A2 <- as.numeric(as.character(data$A2))
    data$A3 <- as.numeric(as.character(data$A3))

    data <- data[row.names(data)!='R2_63',]
    data <- TransQuin(seqs=row.names(data), prop=data$X, A1=data$A1, A2=data$A2, A3=data$A3)

    r <- plotSplice(data, lvl='gene', cname='Ratio', xname='Expected log2 fold change of minor/major', yname='Measured log2 fold change of minor/major')
}

testGenes()