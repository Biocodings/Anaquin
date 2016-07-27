#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

.fitLODR <- function(data,
                     legend='LFC Fold',
                     xlab='Average Counts',
                     ylab='P-value',
                     title='LODR Curves',
                     band='pred',
                     chosenFDR=0.1,
                     multiTest=TRUE,
                     legTitle=NULL,
                     shouldPlotFitting=FALSE)
{
    require(locfit)
    require(qvalue)
    
    stopifnot(!is.null(data$pval))
    stopifnot(!is.null(data$ratio))
    stopifnot(!is.null(data$measured))
    
    data <- data[!is.na(data$pval),]
    data <- data[!is.na(data$measured),]
    
    data$qval <- if(multiTest) qvalue(data$pval)$qvalues else data$pval
    cutoff    <- max(data$pval[data$qval < chosenFDR])
    
    print(paste('FDR threshold:', cutoff))
    
    LODR <- function(pval, mn, cutoff, prob)
    {
        # Eg: 0.02590173 to -1.586671
        cutoff <- log10(cutoff)
        
        # Fit a local regression on the log10 scale
        fit <- locfit(log10(pval)~lp(log10(mn)), maxk=300)
        
        X <- preplot(fit, band=band, newdata=log10(mn))
        
        #plot(fit,band=band,get.data=TRUE,xlim=range(log10(mn)))
        
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
        
        lodr <- segmented.search(fit)
        
        # Bootstrap to estimate uncertainty
        lodr.boot<-NULL
        
        for (ii in 1:2)
        {
            y.boot    <- X$fit + sample(residuals(fit), length(mn))
            fit.boot  <- locfit(y.boot~lp(log10(mn)), maxk=300)
            lodr.boot <- c(lodr.boot, segmented.search(fit.boot)$minimum)
            
            if (ii %% 100 == 0)
            {
                cat ("...")
            }
        }
        
        return (c(lodr$objective, lodr$minimum, quantile(lodr.boot,c(.05,.95))))
    }
    
    lineDat <- NULL;
    
    lodr.resPlot <- NULL; set.seed(1)
    lodr.resLess <- NULL; set.seed(1)
    
    #
    # ---------------------------- Fitting Local Regression ----------------------------
    #
    # Fit a local regression for each ratio. We'll also estimate the confidence interval for each
    # model. Refer to the ERCC paper for more details.
    #
    
    prob <- 0.90
    
    # We'll render for each ratio
    ratios <- sort(data$ratio)
    
    for (ratio in unique(ratios))
    {
        t <- data[data$ratio == ratio,]
        t <- t[t$measured != 0,]
        t <- t[t$pval != 0,]
        
        #
        # 1. Fit a local regression, where each point is approximated by a quadratic function bounded
        #    by a fixed width.
        #
        
        tryCatch (
            {
                print(paste('Estmating LODR for', ratio))
                
                if (shouldPlotFitting)
                {
                    plot(log10(t$measured), log10(t$pval))
                }
                
                # Performs a local regression
                fit <- locfit(log10(t$pval)~lp(log10(t$measured)), maxk=300)
                
                # Plot how the points are fitted
                if (shouldPlotFitting)
                {
                    plot(fit, band=band, get.data=TRUE, main=paste('Local regression for LFC:', ratio))
                }
                
                #
                # Generate new data points and use those data points for the prediction band
                #
                
                x.new <- seq(min(log10(t$measured)), max(log10(t$measured)), length.out=100)
                X <- preplot(fit, band=band, newdata=x.new)
                
                x.new    <- 10^x.new
                fitLine  <- 10^(X$fit)
                fitUpper <- 10^(X$fit + qnorm(prob) * X$se.fit)  # Upper confidence interval
                fitLower <- 10^(X$fit - qnorm(prob) * X$se.fit)  # Lower confidence interval
                
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
                    
                    lodr.resPlot <- rbind(lodr.resPlot, c(round(abs(as.numeric(ratio)), 3), t.res))
                    lodr.resLess <- rbind(lodr.resLess, c(round(abs(as.numeric(ratio)), 3), t.resLess))
                }, error = function(e)
                {
                    print(e)                
                    print(paste('Failed to estimate LODR'))
                })
        }
    }
    
    colnames(lodr.resLess)[1:3] <- c("|log2(Fold)|","MinError","Estimate")
    colnames(lodr.resPlot)[1:3] <- c("ratio","MinError","Estimate")
    colnames(lodr.resLess)[1:3] <- c("ratio","MinError","Estimate")
    
    lodr.resPlot <- as.data.frame(lodr.resPlot)
    lodr.resLess <- as.data.frame(lodr.resLess)
    
    #
    # Construct data for vertical lines for the LODR estimates
    #
    
    arrowDat <- data.frame(ratio = lodr.resPlot$ratio, x = lodr.resPlot[,3], y = cutoff, xend = lodr.resPlot[,3], yend = 0)
    arrowDat$x[grep("<", lodr.resLess[,3])] <- Inf
    arrowDat <- arrowDat[which(is.finite(arrowDat$x)), ]
    
    print('Estimation completed')
    
    lineDat$ratio  <- as.factor(lineDat$ratio)
    arrowDat$ratio <- as.factor(arrowDat$ratio)

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

plotLODR <- function(data,
                     fdr = 0.1,                       
                     title = NULL,
                     xBreaks = NULL,
                     yBreaks = NULL,
                     xLabels = NULL,
                     yLabels = NULL,
                     legTitle = NULL,
                     xlab = NULL,
                     ylab = NULL,
                     shouldFit=TRUE)
{
    data <- data$seqs
    data <- data[data$pval!=0,]
    
    if (shouldFit)
    {
        data <- .fitLODR(data.frame(measured=data$measured, pval=data$pval, ratio=abs(round(data$ratio))))
    }
    else
    {
        data <- data.frame(measured=data$measured, pval=data$pval, ratio=as.factor(abs(round(data$ratio))))
    }

    .plotLODR(data=data,
              xlab=xlab,
              ylab=ylab,
              title=title,
              legTitle=legTitle,
              xBreaks=xBreaks,
              xLabels=xLabels,
              yBreaks=yBreaks,
              yLabels=yLabels)
}