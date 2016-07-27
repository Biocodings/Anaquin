#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

.fitCurve <- function(x, y, algo='locfit', showFitting=TRUE, prob=0.90)
{
    if (showFitting)
    {
        #plot(log10(x), log10(y))
    }

    r <- locfit(y~lp(x), maxk=300)

    if (showFitting)
    {
        plot(r, band='pred', get.data=TRUE, main=paste('Local regression for LFC'))
        #plot(fit,band=band,get.data=TRUE,xlim=range(log10(mn)))
    }

    #
    # Reference: http://www.r-bloggers.com/thats-smooth
    #
    
    x.new <- seq(min(x), max(x), length.out=100)
    X <- preplot(r, band='pred', newdata=x.new)
    x.new <- 10^x.new
    
    uc <- 10^(X$fit + qnorm(prob) * X$se.fit)
    lc <- 10^(X$fit - qnorm(prob) * X$se.fit)

    return (list(fitted=r, uc=uc, lc=lc))
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
    data <- data[!is.na(data$measured),]
    
    if (is.null(data$qval)) { data$qval <- qvalue(data$pval)$qvalues }
    
    # What's the maximum p-value that gives the FDR? This will be the cutoff on the y-axis.
    cutoff <- max(data$pval[data$qval < x$FDR])
    
    LODR <- function(pval, mn, cutoff, prob)
    {
        cutoff <- log10(cutoff)
        
        # Fit a local regression on the log10 scale
        r <- .fitCurve(log10(mn), log10(pval))

        X <- preplot(r$fitted, band=band, newdata=log10(mn))

        find.mn <- function(mn, fit, cutoff, prob)
        {
            X <- preplot(fit, newdata=mn, band=band)
            (X$fit+qnorm(prob)*X$se.fit-cutoff)^2
        }
        
        rng.mn <- range(log10(mn))
        
        # Search in sections to get first crossing
        segmented.search <- function(fit)
        {
            X <- preplot(fit, newdata=min(log10(mn)), band=band)
            
            if ((X$fit + qnorm(prob) * X$se.fit) < cutoff)
            {
                t.lodr <- list(minimum=min(log10(mn)), objective=0)
            }
            else
            {
                ppp<-.2
                t.lodr <- optimize(find.mn, c(rng.mn[1], sum(rng.mn*c(1-ppp,ppp))), fit=fit, cutoff=cutoff, prob=prob)
                
                while (t.lodr$objective > .0001 & ppp<=1)
                {
                    t.lodr<-optimize(find.mn, c(sum(rng.mn*c(1-ppp+.2,ppp-.2)), sum(rng.mn*c(1-ppp,ppp))),
                                     fit=fit, cutoff=cutoff, prob=prob)
                    ppp<-ppp+.2
                }
            }
            
            t.lodr
        }    
        
        lodr <- segmented.search(r$fitted)

        return (c(lodr$objective, lodr$minimum))
    }
    
    lineDat <- NULL;
    
    prob <- 1- x$FDR
    
    # We'll render for each ratio
    ratios <- sort(data$ratio)
    
    for (ratio in unique(ratios))
    {
        t <- data[data$ratio == ratio,]
        t <- t[t$measured != 0,]
        t <- t[t$pval != 0,]
        
        tryCatch (
        {
                print(paste('Estmating LODR for', ratio))
                
                r <- .fitCurve(log10(t$measured),log10(t$pval))
                
                x.new <- seq(min(log10(t$measured)), max(log10(t$measured)), length.out=100)
                X <- preplot(r$fitted, band=band, newdata=x.new)
                
                x.new    <- 10^x.new
                fitLine  <- 10^(X$fit)
                fitUpper <- r$uc
                fitLower <- r$lc

                fitData <- data.frame(x.new, fitLine, fitUpper, fitLower, ratio=ratio)
                lineDat <- rbind(lineDat, fitData)
            }, error = function(e)
            {
                print(e)
                print(paste('Failed to fit a local regression for: ', ratio))
            })
        
        if (ratio != 0)
        {
            tryCatch (
                {
                    t.res <- LODR(t$pval, t$measured, cutoff=cutoff, prob=prob)
                    t.res[-1]<-signif(10^t.res[-1],2)
                    
                    t.resLess <- t.res
                    t.resLess[-1][t.resLess[-1] == signif(min(t$measured),2)] <- paste("<", signif(min(t$measured),2), sep="")
                    t.res[-1][t.res[-1]==signif(min(t$measured),2)] <- Inf
                }, error = function(e)
                {
                    print(e)                
                    print(paste('Failed to estimate LODR'))
                })
        }
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