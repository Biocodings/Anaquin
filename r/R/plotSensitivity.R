#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

.plotSensitivity <- function(data, limit, title='', xlab='', ylab='', showLimit=TRUE, showGuide=FALSE)
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

    if (showLimit)
    {
        r <- min(data[data$y >= limit,]$expected)
        p <- p + geom_vline(xintercept=c(r), linetype="dotted")
        p <- p + geom_label(aes(x=r, y=0.30, label=paste('LOA', r)), colour = "black", show.legend=FALSE)
    }
    
    if (!showGuide)
    {
        p <- p + guides(colour=FALSE)
    }
    
    print(p)
}

plotSensitivity.FusQuin <- function(data, title, xlab, ylab, limit)
{
    .plotSensitivity(data, title=title, xlab=xlab, ylab=ylab, limit=limit)
}

plotSensitivity.TransQuin <- function(data, title, xlab, ylab, limit)
{
    .plotSensitivity(data, title=title, xlab=xlab, ylab=ylab, limit=limit)
}

plotSensitivity <- function(data, title=NULL, xlab=NULL, ylab=NULL, limit=0.98)
{
    stopifnot(class(data) == 'FusQuin' | class(data) == 'TransQuin')

    if (class(data) == 'FusQuin')   { plotSensitivity.FusQuin(data, title=title, xlab=xlab, ylab=ylab, limit=limit) } 
    if (class(data) == 'TransQuin') { plotSensitivity.TransQuin(data, title=title, xlab=xlab, ylab=ylab, limit=limit) } 
}