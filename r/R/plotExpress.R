#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

.plotExpress <- function(data,
                         alpha = 1.0,
                         shouldLog2 = TRUE,
                         showLegend = FALSE,
                         showStats  = "right",
                         title = '',
                         xname = 'Expected log2 fold change of mixture A and B',
                         yname = 'Measured log2 fold change of mixture A and B')
{
    require(ggplot2)

    data     <- data$seqs
    data$x   <- data$expected
    data$y   <- data$measured
    data$grp <- if (shouldLog2) as.factor(abs(log2(data$x))) else as.factor(abs(data$x))
        
    stopifnot(length(data$x) > 0)
    stopifnot(length(data$x) == length((data$y)) || length(data$x) == nrow((data$y)))
    
    if (shouldLog2)
    {
        data$x <- log2(data$x)
        data$y <- log2(data$y)
    }

    isMulti <- is(data$y, 'data.frame')
    
    data$sd   <- if (isMulti) apply(data$y, 1, sd) else NULL
    data$y    <- if (isMulti) rowMeans(data$y)     else data$y
    data$ymax <- if (isMulti) data$y + data$sd     else NULL
    data$ymin <- if (isMulti) data$y - data$sd     else NULL
    
    p <- ggplot(data=data, aes(x=x, y=y)) +
                              xlab(xname) +
                              ylab(yname) +
                           ggtitle(title) +
                geom_point(aes(colour=grp), size=2, alpha=alpha) +
                geom_smooth(method='lm', formula=y ~ x)          +
                labs(colour='Ratio')                             +
                theme_bw()

    y_off <- ifelse(max(data$y) - min(data$y) <= 10, 0.5, 1.0)

    if (showStats == 'right')
    {
        p <- p + annotate("text", label=lm2str(data), x=min(data$x), y=max(data$y)-y_off, size=4, colour='black', parse=TRUE, hjust=1)
    }
    else
    {
        p <- p + annotate("text", label=lm2str(data), x=min(data$x), y=max(data$y)-y_off, size=4, colour='black', parse=TRUE, hjust=0)
    }
    
    p <- p +  theme(axis.title.x=element_text(face='bold', size=12))
    p <- p +  theme(axis.title.y=element_text(face='bold', size=12))

    if (!showLegend)
    {
        p <- p + guides(colour=FALSE)        
    }
    
    if (!is.null(data$sd))
    {
        p <- p + geom_errorbar(aes(ymax=ymax, ymin=ymin), size=0.3, alpha=0.7)
    }
    
    print(p)

    return (list('xname' = xname, 'yname' = yname))
}

plotExpress.MetaQuin <- function(data)
{
    .plotExpress(data, title='Expected expression vs Measured expression',
                 xname='Expected expression (log2)',
                 yname='Measured expression (log2)',
                 showStats='left',
                 showLegend=FALSE)
}

plotExpress.VarQuin <- function(data)
{
    .plotExpress(data, title='Expected expression vs Measured expression',
                       xname='Expected expression (log2)',
                       yname='Measured expression (log2)',
                   showStats='left',
                  showLegend=FALSE)
}

plotExpress.TransQuin <- function(data)
{
    .plotExpress(data, title='Expected expression vs Measured expression',
                       xname='Expected expression (log2)',
                       yname='Measured expression (log2)',
                       showStats='left',
                       showLegend=FALSE)
}

plotExpress <- function(data)
{
    UseMethod("plotExpress", data)
}
