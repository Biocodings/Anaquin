#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

plotFold.TransQuin <- function(data, title)
{
    .plotExpress(data, title=title,
                       showLOQ=FALSE,
                       xlab='Expected log-fold',
                       ylab='Measured log-fold',
                       showStats='left')
}

plotFold <- function(data, title)
{
    stopifnot (class(data) == 'TransQuin')
    
    if (class(data) == 'TransQuin') { plotFold.TransQuin(data, title=title) }
}