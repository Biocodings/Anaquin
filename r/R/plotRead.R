#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

plotRead.TransQuin <- function(data, title, xlab, ylab, showLOQ)
{
    if (is.null(title)) { title <- 'TransQuin Alignment' }
    if (is.null(xlab))  { xlab  <- 'Input concentration (log2 attomol/ul)' }
    if (is.null(ylab))  { ylab  <- 'Reads (log2)' }
    
    .plotExpress(data, title=title, xlab=xlab, ylab=ylab, showLOQ=showLOQ, limitLabel='LOA')
}

plotRead.VarQuin <- function(data, title, xlab, ylab, showLOQ)
{
    if (is.null(title)) { title <- 'VarQuin Alignment' }
    if (is.null(xlab))  { xlab  <- 'Input concentration (log2 attomol/ul)' }
    if (is.null(ylab))  { ylab  <- 'Reads (log2)' }
    
    .plotExpress(data, title=title, xlab=xlab, ylab=ylab, showLOQ=showLOQ, limitLabel='LOA')
}

plotRead <- function(data, title=NULL, xlab=NULL, ylab=NULL, showLOQ=TRUE)
{
    stopifnot(class(data) == 'TransQuin' | class(data) == 'VarQuin')

    if (class(data) == 'TransQuin') { plotRead.TransQuin(data, title=title, xlab=xlab, ylab=ylab, showLOQ=showLOQ) }
    if (class(data) == 'VarQuin')   { plotRead.VarQuin(data, title=title, xlab=xlab, ylab=ylab, showLOQ=showLOQ) }    
}