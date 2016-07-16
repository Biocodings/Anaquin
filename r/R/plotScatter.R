#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

plotScatter <- function(data, xBreaks=NULL, ...)
{
    require(grid)
    require(ggplot2)

    data <- data$seqs

    x <- list(...)
    
    title <- ''
    xlab  <- ''
    ylab  <-
    
    showSD    <- FALSE
    showLOQ   <- TRUE

    if (!is.null(x$xlab))  { xlab  <- x$xlab  }
    if (!is.null(x$ylab))  { ylab  <- x$ylab  }    
    if (!is.null(x$title)) { title <- x$title }
    
    if (!is.null(x$showSD))    { showSD  <- x$showSD      }
    if (!is.null(x$showLOQ))   { showLOQ <- x$showLOQ     }    
    
    if (is.null(x$showInter)) { x$showInter <- FALSE }
        
    # TODO: Fix this....
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
    isMultiSD <- sum(data$sd) > 0 & showSD
    
    isMulti <- isMultiDF | isMultiSD
    
    data$ymax <- NULL
    data$ymin <- NULL

    if (isMultiDF)
    {
        data$sd <- apply(data$y[,-1], 1, sd, na.rm=TRUE)
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
    y_off <- ifelse(max(data$y) - min(data$y) <= 10, 0.7, 0.7)

    if (x$showInter)
    {
        p <- p + geom_vline(xintercept=c(0), linetype='solid', size=0.1)
        p <- p + geom_hline(yintercept=c(0), linetype='solid', size=0.1)
    }
    
    overall <- .lm2str(data)
    above   <- NULL

    loq <- NULL

    if (showLOQ)
    {
        tryCatch({
            # Fit piecewise segmentation
            loq <- showLOQ(data$x, data$y)
        }, error = function(cond)
        {
        })
        
        if (!is.null(loq))
        {
            # Print out the regression above LOQ
            above <- .m2str(loq$model$rModel)
            
            # We can assume the break-point is on the log2-scale. Let's convert it back.
            label <- 2^loq$breaks$k
            
            x <- paste('LOQ:', signif(label, 3))
            x <- paste(x, 'attomol/ul')
            
            p <- p + geom_vline(xintercept=c(loq$breaks$k), linetype='33', size=0.6)
            p <- p + geom_label(aes(x=max(loq$breaks$k), y=min(y)), label=x, colour='black', show.legend=FALSE, hjust=0.1, vjust=0.7)
        }
    }
    
#    title <- list( bquote( paste( "Histogram of " , 'dddd') )  ,
#                   bquote( paste( "Bootstrap samples, Allianz" ) ) )
    
#    t1 <- grid.text(parse(text='bold(Overall)'), gp=gpar(fontsize=11, col="grey24"), draw=TRUE, hjust=1)
#    t2 <- grid.text(parse(text='bold(Above~LOQ)'), gp=gpar(fontsize=11, col="grey24"), draw=TRUE, hjust=1)
    
#    t3 <- grid.text(parse(text=overall), gp=gpar(fontsize=11, col="grey24"), draw=TRUE, hjust=1.0)
#    t4 <- grid.text(parse(text=above), gp=gpar(fontsize=11, col="grey24"), draw=TRUE, hjust=1.0)

#    grid.newpage()
#    g <- arrangeGrob(t1, t3, t2, t4, nrow=2, ncol=2, widths=c(1,1), heights=c(1,1))
#    p <- p + annotation_custom(g, xmin=min(data$x), ymin=max(data$y)-1.5*y_off)
    
    r <- abs(max(data$y) - min(data$y))
    y_off <- 0.06 * r 

    if (showLOQ)
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
    
    if (showLOQ & !is.null(loq))
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

    if (!is.null(xBreaks))
    {
        p <- p + scale_x_continuous(breaks=xBreaks)
    }
    
    if (!is.null(data$sd))
    {
        p <- p + geom_errorbar(aes(ymax=ymax, ymin=ymin), size=0.2, alpha=0.5)
    }
    
    p <- .transformPlot(p)

    #t1 <- grid.text(parse(text=paste(c('bold(Regression)~bold((Overall)): ', overall), collapse='')), draw=TRUE)
    #t2 <- grid.text(parse(text=paste(c('bold(Above~LOQ): ', above), collapse='')), draw=TRUE)

    #grid.newpage()
    #g <- arrangeGrob(t1, t2, nrow=2, heights=c(0.10,0.10))
    #g <- arrangeGrob(t2, nrow=1, heights=c(0.10))    

    #gb <- rectGrob(height = .98, width = .98, gp = gpar(lwd=0, col='black', fill = NA)) # border
    #gt <- gTree(children = gList(g, gb))
    #grid.draw(gt) TODO: Do we need it?
    
#    grid.arrange(p, gt, nrow=2, heights=c(2, 0.22))
#    grid.arrange(p, gt, nrow=2, heights=c(2, 0.15))
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