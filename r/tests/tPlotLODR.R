#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

library('RUnit')
library('Anaquin')

test.PlotLODR_1 <- function()
{
    data('Vignette_5.6.3')
    data <- Vignette_5.6.3

    anaquin <- createAnaquinData(names=row.names(data), measured=data$Mean, ratio=abs(data$ExpLFC), pval=data$Pval)
    
    r <- plotLODR(anaquin, showConf=TRUE)
}

test.PlotLODR_1()