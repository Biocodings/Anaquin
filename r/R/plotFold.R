#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

plotFold.TransQuin <- function(data, x)
{
    .plotExpress(data, title='Expected log-fold vs Measured log-fold',
                       xname='Expected log-fold',
                       yname='Measured log-fold',
                       showStats='left',
                       showLegend=FALSE)
}

plotFold <- function(data, title)
{
    UseMethod("plotFold", data, x=title)
}