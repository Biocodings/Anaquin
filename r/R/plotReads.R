#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

plotReads.TransQuin <- function(data, title, xlab, ylab, showLOQ)
{
    .plotExpress(data, title=title, xlab=xlab, ylab=ylab, showLOQ=showLOQ)
}

plotReads <- function(data, title=NULL, xlab=NULL, ylab=NULL, showLOQ=TRUE)
{
    stopifnot(class(data) == 'TransQuin')

    if (class(data) == 'TransQuin') { plotReads.TransQuin(data, title=title, xlab=xlab, ylab=ylab, showLOQ=showLOQ) }
}