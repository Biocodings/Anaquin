#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

plotScatter <- function(data, ...)
{
    data <- data$seqs

    z <- list(...)
    
    if (is.null(z$xlab))      { z$xlab      <- ''    }
    if (is.null(z$ylab))      { z$ylab      <- ''    }    
    if (is.null(z$title))     { z$title     <- ''    }
    if (is.null(z$showSD))    { z$showSD    <- FALSE }
    if (is.null(z$showLOQ))   { z$showLOQ   <- TRUE  }
    if (is.null(z$showInter)) { z$showInter <- FALSE }

    if (!is.data.frame(data$measured))
    {
        data <- data[!is.na(data$measured),]
        data <- data[!is.infinite(data$measured),]        
    }
    
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
    
    # Should we show standard deviation? Only if we're asked and information available.
    isMultiSD <- sum(data$sd) > 0 & z$showSD

    isMulti <- isMultiDF | isMultiSD
    
    data$ymax <- NULL
    data$ymin <- NULL

    if (isMultiDF)
    {
        data$sd <- apply(data$y, 1, sd, na.rm=TRUE)
        data$y  <- rowMeans(data$y, na.rm=TRUE)
    }
    
    if (isMulti)
    {
        data$ymax <- data$y + data$sd
        data$ymin <- data$y - data$sd
        data <- data[!is.na(data$ymax),]
        data <- data[!is.na(data$ymin),]        
    }
    else
    {
        data$sd <- NULL
    }
    
    data <- data[!is.na(data$y),]
    
    p <- ggplot(data=data, aes_string(x='data$x', y='data$y')) +
                               xlab(z$xlab) +
                               ylab(z$ylab) +
                           ggtitle(z$title) +
                     labs(colour='Ratio') +
                     geom_point(aes_string(colour='grp'), size=2.0, alpha=0.5) +
                geom_smooth(method='lm', formula=y~x, linetype='11', color='black', size=0.5)  +
                theme_bw()

    p <-p + guides(colour=FALSE)
    y_off <- ifelse(max(data$y) - min(data$y) <= 10, 0.7, 0.7)

    if (z$showInter)
    {
        p <- p + geom_vline(xintercept=c(0), linetype='solid', size=0.1)
        p <- p + geom_hline(yintercept=c(0), linetype='solid', size=0.1)
    }
    
    overall <- .lm2str(data)
    above   <- NULL

    LOQ <- NULL

    if (z$showLOQ)
    {
        tryCatch({
            # Fit piecewise segmentation
            LOQ <- showLOQ(data$x, data$y)
        }, error = function(cond)
        {
        })
        
        if (!is.null(LOQ))
        {
            if (LOQ$model$rr > cor(data$x, data$y))
            {
                # Print out the regression above LOQ
                above <- .m2str(LOQ$model$rModel)
                
                # We can assume the break-point is on the log2-scale. Let's convert it back.
                label <- 2^LOQ$breaks$k
                
                t <- paste('LOQ:', signif(label, 3))
                t <- paste(t, 'attomol/ul')
                
                p <- p + geom_vline(xintercept=c(LOQ$breaks$k), linetype='33', size=0.6)
                p <- p + geom_label(aes_string(x='max(LOQ$breaks$k)', y='min(y)'), label=t, colour='black', show.legend=FALSE, hjust=0.1, vjust=0.7)                
            }
            else
            {
                LOQ <- NULL
            }
        }
    }
    
    r <- abs(max(data$y) - min(data$y))
    y_off <- 0.06 * r 

    if (z$showLOQ)
    {
        a <- paste(c('bold(Overall): ', overall), collapse='')
    }
    else
    {
        a <- overall
    }

    overall <- annotate("text",
                        label=a,
                        x=min(data$x),
                        y=max(data$y)-y_off,
                        size=4.0,
                        colour='grey24',
                        parse=TRUE,
                        hjust=0,
                        vjust=0)
    
    p <- p + overall
    
    if (z$showLOQ & !is.null(LOQ))
    {
        above <- annotate("text",
                          label=paste(c('bold(Above)~bold(LOQ): ', above), collapse=''),
                          x=min(data$x),
                          y=max(data$y)-2*y_off,
                          size=4.0,
                          colour='grey24',
                          parse=TRUE,
                          hjust=0,
                          vjust=0)
        p <- p + above
    }

    if (!is.null(z$xBreaks))
    {
        p <- p + scale_x_continuous(breaks=z$xBreaks)
    }
    
    if (!is.null(data$sd))
    {
        p <- p + geom_errorbar(aes_string(ymax='ymax', ymin='ymin'), size=0.2, alpha=0.5)
    }

    p <- .transformPlot(p)
    print(p)

    return (list(LOQ=LOQ$breaks$k))
}