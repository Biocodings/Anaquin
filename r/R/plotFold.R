#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

plotFold.TransQuin <- function(data, title, xlab, ylab, showStats)
{
    if (is.null(title)) { title <- 'TransQuin Differential'  }
    if (is.null(xlab))  { xlab  <- 'Expected log-fold (log2)'}
    if (is.null(ylab))  { ylab  <- 'Measured log-fold (log2)'}    
    
    .plotExpress(data, title=title,
                       showLOQ=FALSE,
                       xlab=xlab,
                       ylab=ylab,
                       showStats=showStats)
}

plotFold.FusQuin <- function(data, title, xlab, ylab, showStats)
{
    if (is.null(title)) { title <- 'VarQuin Differential'    }
    if (is.null(xlab))  { xlab  <- 'Expected log-fold (log2)'}
    if (is.null(ylab))  { ylab  <- 'Measured log-fold (log2)'}    
    
    .plotExpress(data, title=title,
                 showLOQ=FALSE,
                 xlab=xlab,
                 ylab=ylab,
                 showStats=showStats)
}

plotFold <- function(data, title=NULL, xlab=NULL, ylab=NULL, showStats='left')
{
    stopifnot (class(data) == 'TransQuin' || class(data) == 'FusQuin')
    
    if (class(data) == 'TransQuin') { plotFold.TransQuin(data, title=title, xlab=xlab, ylab=ylab, showStats=showStats) }
    if (class(data) == 'FusQuin')   { plotFold.FusQuin(data, title=title, xlab=xlab, ylab=ylab, showStats=showStats)   }    
}