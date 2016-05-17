#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

plotAlign.TransQuin <- function(data, title, xlab, ylab, showLOQ)
{
    if (is.null(title)) { title <- 'TransQuin Alignment' }
    if (is.null(xlab))  { xlab  <- 'Input concentration (log2 attomol/ul)' }
    if (is.null(ylab))  { ylab  <- 'Reads (log2)' }
    
    .plotExpress(data, title=title, xlab=xlab, ylab=ylab, showLOQ=showLOQ, limitLabel='LOA')
}

plotAlign <- function(data, title=NULL, xlab=NULL, ylab=NULL, showLOQ=TRUE)
{
    stopifnot(class(data) == 'TransQuin')

    if (class(data) == 'TransQuin') { plotAlign.TransQuin(data, title=title, xlab=xlab, ylab=ylab, showLOQ=showLOQ) }
}