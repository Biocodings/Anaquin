#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

.getLine <- function(model, y, ratio, band='pred')
{
    x.new <- seq(min(log10(y)), max(log10(y)), length.out=100)
    X <- preplot(model$fitted, band=band, newdata=x.new)
    
    x.new    <- 10^x.new
    fitLine  <- 10^(X$fit)
    fitUpper <- model$uc
    fitLower <- model$lc

    return (data.frame(x.new, fitLine, fitUpper, fitLower, ratio=ratio))
}

LODR <- function(model, pval, mn, cutoff, prob, band)
{
    cutoff <- log10(cutoff)
    
    X <- preplot(model$fitted, band=band, newdata=log10(mn))
    
    lodr <- .segmented.search(model$fitted, mn, band, cutoff, prob, rng.mn=range(log10(mn)))
    
    return (c(lodr$objective, lodr$minimum))
}

.find.mn <- function(mn, fit, cutoff, prob)
{
    X <- preplot(fit, newdata=mn, band='pred')
    (X$fit+qnorm(prob)*X$se.fit-cutoff)^2
}

# Search in sections to get first crossing
.segmented.search <- function(fit, mn, band, cutoff, prob, rng.mn)
{
    X <- preplot(fit, newdata=min(log10(mn)), band=band)
    
    if ((X$fit + qnorm(prob) * X$se.fit) < cutoff)
    {
        t.lodr <- list(minimum=min(log10(mn)), objective=0)
    }
    else
    {
        ppp<-.2
        t.lodr <- optimize(.find.mn, c(rng.mn[1], sum(rng.mn*c(1-ppp,ppp))), fit=fit, cutoff=cutoff, prob=prob)
        
        while (t.lodr$objective > .0001 & ppp<=1)
        {
            t.lodr<-optimize(.find.mn, c(sum(rng.mn*c(1-ppp+.2,ppp-.2)), sum(rng.mn*c(1-ppp,ppp))),
                             fit=fit, cutoff=cutoff, prob=prob)
            ppp<-ppp+.2
        }
    }
    
    t.lodr
}    

.fitCurve <- function(x, y, algo='locfit', showFitting=TRUE, prob=0.90)
{
    if (showFitting)
    {
        #plot(log10(x), log10(y))
    }

    model <- locfit(y~lp(x), maxk=300)

    if (showFitting)
    {
        plot(model, band='pred', get.data=TRUE, main=paste('Local regression for LFC'))
        #plot(fit,band=band,get.data=TRUE,xlim=range(log10(mn)))
    }

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

    band='pred'

    stopifnot(!is.null(data$pval))
    stopifnot(!is.null(data$ratio))
    stopifnot(!is.null(data$measured))
    
    data <- data[!is.na(data$pval),]
    data <- data[data$pval != 0,]    
    data <- data[!is.na(data$measured),]
    
    if (is.null(data$qval)) { data$qval <- qvalue(data$pval)$qvalues }
    
    # What's the maximum p-value that gives the FDR? This will be the cutoff on the y-axis.
    cutoff <- max(data$pval[data$qval < x$FDR])
    
    lineDat <- NULL;
    
    prob <- 1- x$FDR
    
    # We'll render for each ratio
    ratios <- sort(data$ratio)
    
    for (ratio in unique(ratios))
    {
        t <- data[data$ratio == ratio,]

        tryCatch (
        {
            print(paste('Estmating LODR for', ratio))

            # Fitted curve for the ratio
            model <- .fitCurve(log10(t$measured),log10(t$pval))
            
            lineDat <- rbind(lineDat, .getLine(model, t$measured, ratio))
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
        p <- p + geom_line(data=data$lineDat, aes(x=x.new, y=fitLine, colour=ratio), show_guide=FALSE)
        
        if (x$showConf)
        {
            p <- p + geom_ribbon(data=data$lineDat, aes(x=x.new, y=fitLine, ymin=fitLower, ymax=fitUpper, fill=ratio),
                                 alpha=0.3, colour=NA, show_guide=FALSE)
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