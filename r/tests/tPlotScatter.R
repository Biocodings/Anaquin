#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

library('RUnit')
library('Anaquin')

test.PlotScatter_1 <- function()
{
    data(seqCufflinks)
    
    data <- seqCufflinks
    data <- createAnaquinData(names=row.names(data), input=log2(data$InputConcent), measured=log2(data$Observed1))
    
    r <- plotScatter(data, unitTest=TRUE)
    
    checkEqualsNumeric(r$model$breaks, 1.917068888)
}

test.PlotScatter_2 <- function()
{
    data(seqCufflinks)
    
    data <- seqCufflinks
    data <- createAnaquinData(names=row.names(data), input=log2(data$InputConcent), measured=log2(data[,c(2:4)]))
    
    r <- plotScatter(data, unitTest=TRUE)
    
    checkEqualsNumeric(r$model$breaks, 1.917068888)
}

test.PlotScatter_1()
test.PlotScatter_2()