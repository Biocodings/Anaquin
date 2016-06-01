#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

plotJunction.FusQuin <- function(data, title, xlab, ylab)
{
    if (is.null(title)) { title <- 'FusQuin Detection' }
    if (is.null(xlab))  { xlab  <- 'Input concentration (log2 attomol/ul)' }
    if (is.null(ylab))  { ylab  <- 'Junction Reads (log2)' }

    .plotScatter(data, title=title, xlab=xlab, ylab=ylab, showLOQ=FALSE, limitLabel='LOA')
}

plotJunction <- function(data, title=NULL, xlab=NULL, ylab=NULL)
{
    stopifnot(class(data) == 'FusQuin')

    if (class(data) == 'FusQuin') { plotJunction.FusQuin(data, title=title, xlab=xlab, ylab=ylab) }
}