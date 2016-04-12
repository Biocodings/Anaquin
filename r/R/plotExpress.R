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
    data$grp <- abs(data$expected) #ifelse(shouldLog2, abs(log2(data$expect)), abs(data$expect))

    stopifnot(length(data$x) == length((data$y)))
    
    if (length(data$x) == 0)
    {
        warning('Failed to create scatter plot. No data found.')
        return (NULL)
    }
    else if (shouldLog2)
    {
        data$x <- log2(data$x)
        data$y <- log2(data$y)
    }

    data$grp <- as.factor(data$grp)
    
    p <- ggplot(data=data, aes(x=x, y=y)) +
                              xlab(xname) +
                              ylab(yname) +
                           ggtitle(title) +
                geom_point(aes(colour=grp), size=2, alpha=alpha) +
                geom_smooth(method='lm', formula=y ~ x)            +
                labs(colour='Ratio')                               +
                theme_bw()

    x_off <- ifelse(max(data$x) - min(data$x) <= 10, 1.5, 2.0)
    y_off <- ifelse(max(data$y) - min(data$y) <= 10, 0.5, 1.0)

    if (showStats == 'right')
    {
        p <- p + annotate("text", label=lm2str(data), x=min(data$x)+x_off, y=max(data$y)-y_off, size=4, colour='black', parse=TRUE)
    }
    else
    {
        p <- p + annotate("text", label=lm2str(data), x=min(data$x)+x_off, y=max(data$y)-y_off, size=4, colour='black', parse=TRUE)
    }
    
    p <- p +  theme(axis.title.x=element_text(face='bold', size=12))
    p <- p +  theme(axis.title.y=element_text(face='bold', size=12))

    if (!showLegend)
    {
        p <- p + guides(colour=FALSE)        
    }

    print(p)

    return (list('xname' = xname, 'yname' = yname))
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
