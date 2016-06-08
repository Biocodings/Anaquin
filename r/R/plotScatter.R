#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

plotScatter <- function(data, showStats='left', showLOQ=TRUE, title='', xlab='', ylab='', limitLabel='LOQ', xBreaks=NULL)
{
    require(ggplot2)

    data <- data$seqs
    
    # The variable for the x-axis
    data$x <- NULL
    
    # The variable for the y-axis
    data$y <- NULL
    
    if (!is.null(data$input))    { data$x <- data$input    }
    if (!is.null(data$expected)) { data$x <- data$expected }
    if (!is.null(data$measured)) { data$y <- data$measured }
    
    stopifnot(!is.null(data$x))
    stopifnot(!is.null(data$y))
    
    data$grp <- as.factor(abs(data$x))
        
    stopifnot(length(data$x) > 0)
    stopifnot(length(data$x) == length((data$y)) || length(data$x) == nrow((data$y)))
    
    isMulti <- is(data$y, 'data.frame')
    
    data$sd   <- if (isMulti) apply(data$y, 1, sd) else NULL
    data$y    <- if (isMulti) rowMeans(data$y)     else data$y
    data$ymax <- if (isMulti) data$y + data$sd     else NULL
    data$ymin <- if (isMulti) data$y - data$sd     else NULL
    
    p <- ggplot(data=data, aes(x=x, y=y)) +
                               xlab(xlab) +
                               ylab(ylab) +
                           ggtitle(title) +
                geom_point(aes(colour=grp), size=1.0, alpha=1.0) +
                geom_smooth(method='lm', formula=y ~ x)          +
                labs(colour='Ratio')                             +
                theme_bw()
    
    p <-p + guides(colour=FALSE)

    y_off <- ifelse(max(data$y) - min(data$y) <= 10, 0.5, 1.0)

    if (showStats == 'right')
    {
        p <- p + annotate("text", label=lm2str(data), x=min(data$x), y=max(data$y)-y_off, size=4, colour='black', parse=TRUE, hjust=1)
    }
    else
    {
        p <- p + annotate("text", label=lm2str(data), x=min(data$x), y=max(data$y)-y_off, size=4, colour='black', parse=TRUE, hjust=0)
    }
    
    if (showLOQ)
    {
        loq <- showLOQ(data$x, data$y)
        p <- p + geom_vline(xintercept=c(loq$breaks$k), linetype="dotted")
        p <- p + geom_label(aes(x=loq$breaks$k, y=0, label=paste(limitLabel, signif(loq$breaks$k, 3))), colour = "black", show.legend=FALSE)
    }

    if (!is.null(xBreaks))
    {
        p <- p + scale_x_continuous(breaks=xBreaks)
    }
    
    p <- p + theme(axis.title.x=element_text(face='bold', size=12))
    p <- p + theme(axis.title.y=element_text(face='bold', size=12))

    if (!is.null(data$sd))
    {
        p <- p + geom_errorbar(aes(ymax=ymax, ymin=ymin), size=0.3, alpha=0.7)
    }

    print(p)
}

plotExpress.TransQuin <- function(data, title, xlab, ylab, showLOQ)
{
    # TODO: Fix this
    if (is.null(title)) { title <- 'Isoform Expression' }
    if (is.null(xlab))  { xlab <- 'Input concentration (log2) '}
    if (is.null(ylab))  { ylab <- 'FPKM (log2) '}
    
    xBreaks <- c(-3, 0, 6, 9, 12, 15)
    .plotExpress(data, title=title, xlab=xlab, ylab=ylab, xBreaks=xBreaks, showLOQ=showLOQ)
}