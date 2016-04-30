#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

plotFold.TransQuin <- function(data, title, xlab, ylab, showStats)
{
    .plotExpress(data, title=title,
                       showLOQ=FALSE,
                       xlab=xlab,
                       ylab=ylab,
                       showStats=showStats)
}

plotFold.FusQuin <- function(data, title, xlab, ylab, showStats)
{
    .plotExpress(data, title=title,
                 showLOQ=FALSE,
                 xlab=xlab,
                 ylab=ylab,
                 showStats=showStats)
}

plotFold <- function(data, title='', xlab='', ylab='', showStats='left')
{
    stopifnot (class(data) == 'TransQuin' || class(data) == 'FusQuin')
    
    if (class(data) == 'TransQuin') { plotFold.TransQuin(data, title=title, xlab=xlab, ylab=ylab, showStats=showStats) }
    if (class(data) == 'FusQuin')   { plotFold.FusQuin(data, title=title, xlab=xlab, ylab=ylab, showStats=showStats)   }    
}