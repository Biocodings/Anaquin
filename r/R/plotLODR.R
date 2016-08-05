#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

.getLODR <- function(ratio, model, x, y, pval)
{
    knots <- seq(min(x), max(x), length.out=5000)

    preds <- predict(model, band='pred', newdata=knots)
    preds <- 10^preds$fit
    
    # Does the curve intersect the y-axis?
    if (min(preds) > pval | ratio == 0)
    {
        return (Inf)
    }
    
    # Difference between predictions and p-value cutoff (y-axis)
    diff <- abs(preds - pval)

    #
    # LODR is simply the point that is closest to the intersection. This is not exact solution
    # but LODR itself is a rough estimate. It's not possible for root solving because local
    # regression is non-parametric. The ERCC package uses a complicated estimation algorithm,
    # but there is really no need for that.
    #

    return (10^knots[which.min(diff)])
}

.fitCurve <- function(ratio, x, y, pval, algo='locfit', showFitting=FALSE)
{
    if (algo == 'locfit')
    {
        model <- locfit(y~lp(x))

        # Points where the curve is approximated
        knots <- seq(min(x), max(x), length.out=100)
        
        # Predictions for the knots
        kpred <- predict(model, band='pred', newdata=knots)

        uc <- 10^(kpred$fit + qnorm(pval) * kpred$se.fit)
        lc <- 10^(kpred$fit - qnorm(pval) * kpred$se.fit)
    }

    if (showFitting) { plot(model, band='pred', get.data=TRUE) }

    data <- data.frame(ratio=ratio, knots=10^knots, pred=10^kpred$fit, uc=uc, lc=lc)

    if (any(!is.finite(data$lc)) | any(!is.finite(data$uc)))
    {
        data <- data.frame()
    }

    return (list(LODR=.getLODR(ratio, model, x, y, pval), data=data))
}

.fitLODR <- function(data, ...)
{
    x <- list(...)

    stopifnot(!is.null(data$pval))
    stopifnot(!is.null(data$ratio))
    stopifnot(!is.null(data$measured))
    
    data <- data[!is.na(data$pval),]
    data <- data[data$pval != 0,]    
    data <- data[!is.na(data$measured),]
    
    if (is.null(x$FDR))     { x$FDR <- 0.05 }
    if (is.null(data$qval)) { data$qval <- qvalue(data$pval, fdr.level=x$FDR)$qvalues }

    data <- data[!is.na(data$qval),]
    
    # What's the maximum p-value that gives the FDR? This will be the cutoff on the y-axis.
    pval <- max(data$pval[data$qval < x$FDR])

    print(paste(c('Threshold: ', pval), collapse=''))

    curve <- NULL
    limit <- NULL
    
    for (ratio in unique(sort(data$ratio)))
    {
        t <- data[data$ratio == ratio,]

        tryCatch (
        {
            print(paste('Estmating LODR for', ratio))
            
            r <- .fitCurve(ratio=ratio, log10(t$measured), log10(t$pval), pval=pval)
            
            curve <- rbind(curve, r$data)
            limit <- rbind(limit, data.frame(Ratio=ratio, LODR=r$LODR))
        }, error = function(e)
        {
            print(e)
            print(paste('Failed to curve fit for: ', ratio))
        })
    }
    
    curve$ratio <- as.factor(curve$ratio)

    return(list(measured=data$measured,
                pval=data$pval,
                ratio=as.factor(data$ratio),
                curve=curve,
                limit=limit))
}

.plotLODR <- function(data, ...)
{
    stopifnot(!is.null(data$pval))
    stopifnot(!is.null(data$ratio))
    stopifnot(!is.null(data$measured))
    
    x <- list(...)
    
    if (is.null(x$size))     { x$size     <- 3       }
    if (is.null(x$legTitle)) { x$legTitle <- 'Ratio' }
    if (is.null(x$showConf)) { x$showConf <- FALSE   }

    df <- data.frame(measured=data$measured, pval=data$pval, ratio=data$ratio)
    p  <- ggplot(df, aes_string(x='measured', y='pval', colour='ratio'))   +
                         geom_point(size=x$size, alpha=0.5)   +
                         labs(colour=x$legTitle)              +
                         theme_bw()

    if (!is.null(x$xlab))     { p <- p + xlab(x$xlab)     }
    if (!is.null(x$ylab))     { p <- p + ylab(x$ylab)     }
    if (!is.null(x$title))    { p <- p + ggtitle(x$title) }
    if (!is.null(x$legTitle)) { p <- p + labs(colour=x$legTitle) }

    p <- p + geom_line(data=data$curve, aes_string(x='knots', y='pred', colour='ratio'), show.legend=FALSE)

    if (x$showConf & !is.null(data$curve))
    {
        p <- p + geom_ribbon(data=data$curve, aes_string(x='knots', y='pred', ymin='lc', ymax='uc',
                             fill='ratio'), alpha=0.3, colour=NA, show.legend=FALSE)
    }

    if (!is.null(x$xBreaks))
    {
        p <- p + scale_x_log10(breaks=x$xBreaks)
    }
    else
    {
        p <- p + scale_x_log10()
    }

    if (!is.null(x$yBreaks))
    {
        p <- p + scale_y_log10(breaks=x$yBreaks)
    }
    else
    {
        p <- p + scale_y_log10()        
    }
    
    if (!is.null(data$limit))
    {
        limit <- data$limit[is.finite(data$limit$LODR),]
        print(kable(limit))

        limit <- cbind(limit, y=c(0))
        limit <- cbind(limit, yend=c(1))        
        
        p <- p + geom_segment(data=limit, aes_string(x='LODR', y='y', xend='LODR', yend='yend', color='as.factor(limit$Ratio)'), linetype='33', size=0.6)
    }
    
    print(.transformPlot(p))
}

plotLODR <- function(data, title, xlab, ylab, FDR, legTitle, showConf)
{
    data <- data$seqs
    data <- data[data$pval!=0,]
    
    x <- list(...)
    
    # Sequin ratio groups
    ratio <- as.factor(abs(round(data$ratio)))
    
    # Measured variable (eg: average counts)
    measured <- data$measured
    
    # Measured p-value probability
    pval <- data$pval

    if (!is.null(data$qval))
    {
        data <- .fitLODR(data.frame(measured=data$measured,
                                    title=x$title,
                                    xlab=x$xlab,
                                    ylab=x$ylab,
                                    FDR=x$FDR,
                                    legTitle=x$legTitle,
                                    showConf=x$showConf,
                                    pval=data$pval,
                                    qval=data$qval,
                                    ratio=abs(round(data$ratio))), ...)
    }
    else
    {
        data <- .fitLODR(data.frame(measured=data$measured,
                                    title=x$title,
                                    xlab=x$xlab,
                                    ylab=x$ylab,
                                    FDR=x$FDR,
                                    legTitle=x$legTitle,
                                    showConf=x$showConf,
                                    pval=data$pval,
                                    ratio=abs(round(data$ratio))), ...)
    }

    .plotLODR(data=data, ...)
}