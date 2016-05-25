#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

plotKExpress.MetaQuin <- function(data, title, xlab, ylab, showLOQ)
{
    title <- ifelse(is.null(title), 'Metagenomics Detection', title)
    xlab  <- ifelse(is.null(xlab),  'Input Concentration (log2 attomoles/ul)', xlab)
    ylab  <- ifelse(is.null(ylab),  'K-mer Coverage (log2)', ylab)

    .plotScatter(data, title=title, xlab=xlab, ylab=ylab, showLOQ=showLOQ)
}

plotKExpress.VarQuin <- function(data, title, xlab, ylab, showLOQ)
{
    title <- ifelse(is.null(title), 'VarQuin Detection', title)
    xlab  <- ifelse(is.null(xlab),  'Input Concentration (log2 attomoles/ul)', xlab)
    ylab  <- ifelse(is.null(ylab),  'Gene Abundance FPKM (log2)', ylab)
    
    .plotScatter(data, title=title, xlab=xlab, ylab=ylab, showLOQ=showLOQ)
}

plotKExpress <- function(data, title=NULL, xlab=NULL, ylab=NULL, showLOQ=TRUE)
{
    stopifnot(class(data) == 'MetaQuin' | class(data) == 'TransQuin' | class(data) == 'VarQuin')

    if (class(data) == 'MetaQuin') { plotKExpress.MetaQuin(data, title=title, xlab=xlab, ylab=ylab, showLOQ=showLOQ)  }
    if (class(data) == 'VarQuin')  { plotKExpress.VarQuin(data, title=title, xlab=xlab, ylab=ylab, showLOQ=showLOQ)  }    
}