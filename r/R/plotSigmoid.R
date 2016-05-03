#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

.plotSigmoid <- function(data, title='', xlab='', ylab='', showGuide=FALSE)
{
    require(ggplot2)

    data     <- data$seqs
    data$x   <- data$expected
    data$y   <- data$measured
    data$grp <- as.factor(round(abs(data$x)))

    stopifnot(length(data$x) > 0)
    stopifnot(length(data$x) == length((data$y)) || length(data$x) == nrow((data$y)))
    
    p <- ggplot(data=data, aes(x=x, y=y)) +
                               xlab(xlab) +
                               ylab(ylab) +
                               ggtitle(title) +
                               geom_point(aes(colour=grp), size=2, alpha=1.0) +        
                               theme_bw()

    p <- p + theme(axis.title.x=element_text(face='bold', size=12))
    p <- p + theme(axis.title.y=element_text(face='bold', size=12))

    if (!showGuide)
    {
        p <- p + guides(colour=FALSE)
    }
    
    print(p)
}

plotSigmoid.FusQuin <- function(data, title, xlab, ylab)
{
    .plotSigmoid(data, title=title, xlab=xlab, ylab=ylab)
}

plotSigmoid <- function(data, title=NULL, xlab=NULL, ylab=NULL, showLOQ=TRUE)
{
    stopifnot(class(data) == 'FusQuin')

    if (class(data) == 'FusQuin') { plotSigmoid.FusQuin(data, title=title, xlab=xlab, ylab=ylab) } 
}