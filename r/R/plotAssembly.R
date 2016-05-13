#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

plotAssembly.TransQuin <- function(data, title, xlab, ylab, limit)
{
    .plotSensitivity(data, title=title, xlab=xlab, ylab=ylab, limit=limit)
}

plotAssembly <- function(data,
                         title='Assembly Sensitivity',
                         xlab='Input Concentration (log2 attomole/ul)',
                         ylab='Sensitivity',
                         limit=0.98)
{
    stopifnot(class(data) == 'TransQuin')

    if (class(data) == 'TransQuin') { plotAssembly.TransQuin(data, title=title, xlab=xlab, ylab=ylab, limit=limit) } 
}