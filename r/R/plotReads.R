#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

plotReads.TransQuin <- function(data, title, xlab, ylab, showLOQ)
{
    if (is.null(title)) { title <- 'TransQuin Reads' }
    if (is.null(xlab))  { xlab  <- 'Input concentration (log2 attomol/ul)' }
    if (is.null(ylab))  { ylab  <- 'Reads (log2)' }
    
    .plotExpress(data, title=title, xlab=xlab, ylab=ylab, showLOQ=showLOQ)
}

plotReads <- function(data, title=NULL, xlab=NULL, ylab=NULL, showLOQ=TRUE)
{
    stopifnot(class(data) == 'TransQuin')

    if (class(data) == 'TransQuin') { plotReads.TransQuin(data, title=title, xlab=xlab, ylab=ylab, showLOQ=showLOQ) }
}