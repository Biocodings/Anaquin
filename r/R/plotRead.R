#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

plotRead.MetaQuin <- function(data, title, xlab, ylab, showLOQ)
{
    title <- ifelse(is.null(title), 'Metagenomics Detection', title)
    xlab  <- ifelse(is.null(xlab),  'Input Concentration (log2 attomoles/ul)', xlab)
    ylab  <- ifelse(is.null(ylab),  'FPKM (log2)', ylab)
    
    .plotScatter(data, title=title, xlab=xlab, ylab=ylab, showLOQ=showLOQ)
}

plotRead.TransQuin <- function(data, title, xlab, ylab, showLOQ)
{
    if (is.null(title)) { title <- 'TransQuin Alignment' }
    if (is.null(xlab))  { xlab  <- 'Input concentration (log2 attomol/ul)' }
    if (is.null(ylab))  { ylab  <- 'Reads (log2)' }

    .plotScatter(data, title=title, xlab=xlab, ylab=ylab, showLOQ=showLOQ, limitLabel='LOA')
}

plotRead <- function(data, title=NULL, xlab=NULL, ylab=NULL, showLOQ=TRUE)
{
    stopifnot(class(data) == 'TransQuin' | class(data) == 'MetaQuin')
    
    # Data required for number of reads aligned for each sequin
    stopifnot(!is.null(data$seqs$measured))
    
    if (class(data) == 'TransQuin') { plotRead.TransQuin(data, title=title, xlab=xlab, ylab=ylab, showLOQ=showLOQ) }
    if (class(data) == 'MetaQuin')  { plotRead.MetaQuin(data, title=title, xlab=xlab, ylab=ylab, showLOQ=showLOQ) }    
}