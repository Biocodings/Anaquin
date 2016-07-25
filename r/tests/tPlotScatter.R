#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

library('RUnit')
library('Anaquin')

test.PlotScatter_1 <- function()
{
    data("Vignette_5.4.6.1")
    
    data <- Vignette_5.4.6.1
    data <- CreateDataForAnaquin(names=row.names(data), input=log2(data$InputConcent), measured=log2(data$Observed))
    
    r <- plotScatter(data)
    
    checkEqualsNumeric(r$LOQ, 1.917068888)
}

test.PlotScatter_1()