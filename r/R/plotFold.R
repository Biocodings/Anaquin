#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

plotLogFold.TransQuin <- function(data)
{
    .plotExpress(data, title='Expected log-fold vs Measured log-fold',
                       xname='Expected log-fold',
                       yname='Measured log-fold',
                       shouldLog2=FALSE,
                       showStats='left',
                       showLegend=FALSE)
}

plotLogFold <- function(data)
{
    UseMethod("plotLogFold", data)
}
