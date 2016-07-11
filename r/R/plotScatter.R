#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

plotScatter <- function(data, showIntercept=FALSE, showLOQ=TRUE, title='', xlab='', ylab='', xBreaks=NULL)
{
    require(ggplot2)

    data <- data$seqs
    
    # The variable for the x-axis
    data$x <- NULL
    
    # The variable for the y-axis
    data$y <- NULL
    
    if (!is.null(data$input))
    { 
        data$expected <- data$input
    }

    if (!is.null(data$expected)) { data$x <- data$expected }
    if (!is.null(data$measured)) { data$y <- data$measured }
    
    stopifnot(!is.null(data$x))
    stopifnot(!is.null(data$y))
    
    data <- data[!is.na(data$expected),]
    data$grp <- as.factor(abs(data$x))
        
    stopifnot(length(data$x) > 0)
    stopifnot(length(data$x) == length((data$y)) || length(data$x) == nrow((data$y)))
    
    isMultiDF <- is(data$y, 'data.frame')
    isMultiSD <- sum(data$sd) > 0
    isMulti   <- isMultiDF | isMultiSD
    
    data$ymax <- NULL
    data$ymin <- NULL

    if (isMultiDF)
    {
        data$sd <- apply(data$y, 1, sd)
        data$y  <- rowMeans(data$y)
    }
    
    if (isMulti)
    {
        data$ymax <- data$y + data$sd
        data$ymin <- data$y - data$sd        
    }
    else
    {
        data$sd <- NULL
    }
    
    #data$sd   <- if (isMulti) apply(data$y, 1, sd) else NULL
    #data$y    <- if (isMulti) rowMeans(data$y)     else data$y
    #data$ymax <- if (isMulti) data$y + data$sd     else NULL
    #data$ymin <- if (isMulti) data$y - data$sd     else NULL

    data <- data[!is.na(data$y),]
    
    p <- ggplot(data=data, aes(x=x, y=y)) +
                               xlab(xlab) +
                               ylab(ylab) +
                           ggtitle(title) +
                     labs(colour='Ratio') +
                     geom_point(aes(colour=grp), size=2.0, alpha=0.5) +
                geom_smooth(method='lm', formula=y~x, linetype='11', color='black', size=0.5)  +
                theme_bw()

    p <-p + guides(colour=FALSE)
    y_off <- ifelse(max(data$y) - min(data$y) <= 10, 0.5, 0.5)

    # Position the regression model to the left
    p <- p + annotate("text",
                      label=lm2str(data),
                      x=min(data$x),
                      y=max(data$y)-y_off, size=4, colour='black', parse=TRUE, hjust=0)

    if (showIntercept)
    {
        p <- p + geom_vline(xintercept=c(0), linetype='solid', size=0.1)
        p <- p + geom_hline(yintercept=c(0), linetype='solid', size=0.1)
    }

    if (showLOQ)
    {
        # Fit piecewise segmentation
        loq <- showLOQ(data$x, data$y)

        # We can assume the break-point is on the log2-scale. Let's convert it back.
        label <- 2^loq$breaks$k
        
        x <- paste('LOQ:', signif(label, 3))
        x <- paste(x, 'attomol/ul')

        p <- p + geom_vline(xintercept=c(loq$breaks$k), linetype='33', size=0.6)
        p <- p + geom_label(aes(x=loq$breaks$k, y=min(y)), label=x, colour='black', show.legend=FALSE)
    }

    if (!is.null(xBreaks))
    {
        p <- p + scale_x_continuous(breaks=xBreaks)
    }
    
    if (!is.null(data$sd))
    {
        p <- p + geom_errorbar(aes(ymax=ymax, ymin=ymin), size=0.3, alpha=0.7)
    }
    
    p <- .transformPlot(p)
    print(p)
}

plotFold <- function(data, showIntercept=FALSE, title='', xlab='', ylab='', xBreaks=NULL)
{
    plotScatter(data, showIntercept=showIntercept, showLOQ=FALSE, title=title, xlab=xlab, ylab=ylab, xBreaks=xBreaks)
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