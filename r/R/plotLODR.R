#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

.smoothCurve <- function(model, y, ratio, algo='locfit', band='pred')
{
    if (algo == 'locfit')
    {
        x.new <- seq(min(log10(y)), max(log10(y)), length.out=100)

        x <- predict(model$fitted, newdata=x.new)
        
        x.new    <- 10^x.new
        fitLine  <- 10^(x)
        fitUpper <- model$uc
        fitLower <- model$lc
    }
    else if (algo == 'loess')
    {
        
    }

    return (data.frame(x.new, fitLine, fitUpper, fitLower, ratio=ratio))
}

.getLODR <- function()
{
    
}

.fitCurve <- function(x, y, algo='locfit', showFitting=FALSE, prob=0.90)
{
    if (algo == 'locfit')
    {
        model <- locfit(y~lp(x), maxk=300)
    }
    else if (algo == 'loess')
    {
        
    }

    if (showFitting) { plot(model, band='pred', get.data=TRUE) }

    #
    # Reference: http://www.r-bloggers.com/thats-smooth
    #
    
    x.new <- seq(min(x), max(x), length.out=100)
    X <- preplot(model, band='pred', newdata=x.new)
    x.new <- 10^x.new
    
    uc <- 10^(X$fit + qnorm(prob) * X$se.fit)
    lc <- 10^(X$fit - qnorm(prob) * X$se.fit)

    return (list(fitted=model, uc=uc, lc=lc))
}

.fitLODR <- function(data, ...)
{
    require(locfit)
    require(qvalue)
    
    x <- list(...)

    stopifnot(!is.null(data$pval))
    stopifnot(!is.null(data$ratio))
    stopifnot(!is.null(data$measured))
    
    data <- data[!is.na(data$pval),]
    data <- data[data$pval != 0,]    
    data <- data[!is.na(data$measured),]
    
    if (is.null(data$qval)) { data$qval <- qvalue(data$pval)$qvalues }
    
    # What's the maximum p-value that gives the FDR? This will be the cutoff on the y-axis.
    cutoff <- max(data$pval[data$qval < x$FDR])

    # We'll render for each ratio
    ratios <- sort(data$ratio)

    lineDat <- NULL;

    for (ratio in unique(ratios))
    {
        t <- data[data$ratio == ratio,]

        tryCatch (
        {
            print(paste('Estmating LODR for', ratio))

            # Fitted curve for the ratio
            model <- .fitCurve(log10(t$measured),log10(t$pval))
            
            lineDat <- rbind(lineDat, .smoothCurve(model, t$measured, ratio))
        }, error = function(e)
        {
            print(e)
            print(paste('Failed to curve fit for: ', ratio))
        })
    }
    
    lineDat$ratio <- as.factor(lineDat$ratio)

    return(list(measured=data$measured,
                pval=data$pval,
                ratio=as.factor(data$ratio),
                lineDat=lineDat))
}

.plotLODR <- function(data, ...)
{
    require(qvalue)
    require(ggplot2)

    stopifnot(!is.null(data$pval))
    stopifnot(!is.null(data$ratio))
    stopifnot(!is.null(data$measured))
    
    x <- list(...)
    
    if (is.null(x$size))     { x$size     <- 3       }
    if (is.null(x$showConf)) { x$showConf <- TRUE    }
    if (is.null(x$legTitle)) { x$legTitle <- 'Ratio' }
    
    df <- data.frame(measured=data$measured, pval=data$pval, ratio=data$ratio)
    
    p <- ggplot(df, aes(x=measured, y=pval, colour=ratio)) +
                      geom_point(size=x$size, alpha=0.5)   +
                      labs(colour=x$legTitle)              +
                      theme_bw()

    if (is.null(x$xBreaks))
    {
        x$xBreaks <- 10^(0:round(log10(max(df$measured)))-1)
    }
    
    if (is.null(x$xLabels))
    {
        x$xLabels <- paste('1e+0', 0:(length(x$xBreaks)-1), sep='')
    }
    
    if (is.null(x$yBreaks))
    {
        #x$yBreaks <- c(min(df$pval), 1e-300, 1e-200, 1e-100, 1e-10, 1)
        x$yBreaks <- c(min(df$pval), 1e-100, 1e-80, 1e-60, 1e-40, 1e-20, 1.00)
    }
    
    if (!is.null(x$xlab))     { p <- p + xlab(x$xlab)     }
    if (!is.null(x$ylab))     { p <- p + ylab(x$ylab)     }
    if (!is.null(x$title))    { p <- p + ggtitle(x$title) }
    if (!is.null(x$legTitle)) { p <- p + labs(colour=x$legTitle) }
    
    if (!is.null(data$lineDat))
    {
        p <- p + geom_line(data=data$lineDat, aes(x=x.new, y=fitLine, colour=ratio), show.legend=FALSE)
        
        if (x$showConf)
        {
            p <- p + geom_ribbon(data=data$lineDat, aes(x=x.new, y=fitLine, ymin=fitLower, ymax=fitUpper, fill=ratio),
                                 alpha=0.3, colour=NA, show.legend=FALSE)
        }
    }
    
    if (!is.null(x$arrowDat))
    {
        p <- p + geom_segment(data=x$arrowDat, aes(x=x, y=y, xend=xend, yend=yend, colour=ratio), 
                              lineend="round", arrow=grid::arrow(length = grid::unit(0.5, 'cm')), size=2, alpha=0.6) +
            labs(colour='Log-Fold') +
            p <- p + geom_hline(yintercept=x$cutoff, linetype=2, size=2)
    }
    
    if (!is.null(x$xMin) & !is.null(x$xMax))
    {
        p <- p + xlim(x$xMin, x$xMax)
    }
    
    if (!is.null(x$yMin) & !is.null(x$yMax))
    {
        p <- p + ylim(x$yMin, x$yMax)
    }
    
    if (!is.null(x$xBreaks) & !is.null(x$xLabels))
    {
        p <- p + scale_x_log10(breaks=x$xBreaks, labels=x$xLabels)
    }
    else
    {
        p <- p + scale_x_log10(limits=c(min(data$measured), max(data$measured)))
        p <- p + scale_x_log10(limits=c(min(data$measured), max(data$measured)), breaks=c(arrowDat$x, round(max(data$baseMean))))
    }
    
    if (!is.null(x$yBreaks))
    {
        p <- p + scale_y_log10(breaks=x$yBreaks)
    }
    
    p <- .transformPlot(p)
    print(p)
}

plotLODR <- function(data, ...)
{
    data <- data$seqs
    data <- data[data$pval!=0,]
    
    # Should we attempt for curve fitting? (not every distribution can be fitted)
    shouldFit <- TRUE

    x <- list(...)
    
    if (!is.null(x$shouldFit)) { shouldFit <- x$shouldFit }
    
    # Sequin ratio groups
    ratio <- as.factor(abs(round(data$ratio)))
    
    # Measured variable (eg: average counts)
    measured <- data$measured
    
    # Measured p-value probability
    pval <- data$pval

    if (shouldFit)
    {
        data <- .fitLODR(data.frame(measured=data$measured,
                                    pval=data$pval,
                                    ratio=abs(round(data$ratio))), ...)
    }
    else
    {
        data <- data.frame(measured=data$measured,
                           pval=data$pval,
                           ratio=as.factor(abs(round(data$ratio))))
    }

    .plotLODR(data=data, ...)
}