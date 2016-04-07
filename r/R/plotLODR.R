#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

.plotLODR <- function(data,
                      band='pred',
                      chosenFDR=0.1,
                      shouldTable=FALSE,
                      multiTest=TRUE)
{
    require(locfit)

    stopifnot(!is.null(data$ratio) & !is.null(data$pval) & !is.null(data$abund))
    
    data <- data[!is.na(data$pval),]
    data <- data[!is.na(data$abund),]

    #
    # Estimate the q-value for false discovery rate.
    #
    #    https://github.com/Bioconductor-mirror/erccdashboard/blob/631b691a51db54cb165ff2de1f9a5e06608347bd/R/geneExprTest.R
    #
    
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
    
    # Apply to loaded data
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
        t <- t[t$abund != 0,]
        t <- t[t$pval != 0,]

        #
        # 1. Fit a local regression, where each point is approximated by a quadratic function bounded
        #    by a fixed width.
        #
        
        tryCatch (
        {
            print(paste('Estmating LODR for', ratio))
            
            plot(log10(t$abund), log10(t$pval))
            
            # Performs a local regression
            fit <- locfit(log10(t$pval)~lp(log10(t$abund)), maxk=300)
            
            # Plot how the points are fitted
            plot(fit, band=band, get.data=TRUE, main=paste('Local regression for LFC:', ratio))
            
            #
            # Generate new data points and use those data points for the prediction band
            #
            
            x.new <- seq(min(log10(t$abund)), max(log10(t$abund)), length.out=100)
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
                t.res <- LODR(t$pval, t$abund, cutoff=cutoff, prob=prob)
                t.res[-1]<-signif(10^t.res[-1],2)
                
                #if (t.res[1]>.01)
                #{
                #   t.res[2]<-Inf
                #  t.res[3:4]<-NA
                #}
                
                t.resLess <- t.res
                t.resLess[-1][t.resLess[-1] == signif(min(t$abund),2)] <- paste("<", signif(min(t$abund),2), sep="")
                t.res[-1][t.res[-1]==signif(min(t$abund),2)] <- Inf

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
    
    if (shouldTable)
    {
        legendLabels <- c('4', '3', '2', '1') # TODO: Fix this
        
        lodr.resPlot$ratio <- as.character(legendLabels)
        lodr.resLess$ratio <- as.character(legendLabels) 
        annoTable <- lodr.resLess[-c(2,4,5)]
        
        colnames(annoTable) <- c("Ratio", expression("LODR Estimate"))
        cat("\n")
        print(annoTable, quote = FALSE, rofw.names = FALSE)    
        
        my_table <- tableGrob(d=annoTable, rows=NULL)
    }
    
    # What color scheme should we use?
    #cols <- colors(length(data$ratio))
    
    #data$eLogLF    <- as.factor(data$ratio)
    lineDat$ratio  <- as.factor(lineDat$ratio)
    arrowDat$ratio <- as.factor(arrowDat$ratio)
    
    x <- data.frame(abund=data$abund, pval=data$pval, ratio=as.factor(data$ratio))

    plotLODR.Plot(data=x, lineDat=lineDat, shouldBand=TRUE, cutoff=cutoff)
}

plotAllelePValue <- function(data, ..., xBreaks=c(-3, -2, -1, 0))
{
    require(plyr)
    
    data$ratio <- revalue(data$ratio, c('1'='FP'))
    
    xLabels <- xBreaks
    xLabels[xLabels==0] <- 'FP'
    
    plotLODR.Plot(data, title='Expected allele frequency vs p-value (SNP)', xname='Expected allele frequency (log10)', yname='P-value (log10)', xBreaks=xBreaks, xLabels=xLabels)
}

#
# Construct the LODR plot. The x-axis would be the measured expression while the y-axis would be the p-value.
# This function doesn't apply any transformation to the inputs.
#

plotLODR.Plot <- function(data, ...)
{
    require(grid)
    require(qvalue)
    require(ggplot2)
    require(gridExtra)
    
    x <- list(...)
    p <- ggplot(data, aes(x=abund, y=pval, colour=ratio)) + geom_point(size=3) + theme_bw()

    
    x$xBreaks <- c(1, 10, 100, 1000, 10000)
    x$xLabels <- c('1e+00', '1e+01', '1e+02', '1e+03', '1e+04')
    x$yBreaks <- c(1e-300, 1e-200, 1e-100, 1e-10, 1.00)
    
    
    if (!is.null(x$xname)) { p <- p + xlab(x$xname)    }
    if (!is.null(x$yname)) { p <- p + ylab(x$yname)    }
    if (!is.null(x$title)) { p <- p + ggtitle(x$title) }

    if (!is.null(x$lineDat))
    {
        p <- p + geom_line(data=x$lineDat, aes(x=x.new, y=fitLine, colour=ratio), show_guide=FALSE)
        
        if (x$shouldBand)
        {
            p <- p + geom_ribbon(data=x$lineDat, aes(x=x.new, y=fitLine, ymin=fitLower, ymax=fitUpper, fill=ratio),
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
        #p <- p + scale_x_continuous(breaks=x$xBreaks, labels=x$xLabels)        
    }
    else
    {
        p <- p + scale_x_log10(limits=c(min(data$abund), max(data$abund)))
        p <- p + scale_x_log10(limits=c(min(data$abund), max(data$abund)), breaks=c(arrowDat$x, round(max(data$baseMean))))
    }
    
    if (!is.null(x$yBreaks))
    {
        p <- p + scale_y_log10(breaks=x$yBreaks)
    }
    
    if (!is.null(x$drawTable))
    {
        annotLODRplot <- grid.arrange(arrangeGrob(grobs = list(p, my_table), ncol = 1, heights = c(2,0.5)))
    }
    
    print(p)
}

plotLODR.VarQuin <- function(data, title=NULL)
{
    # Probability under the null hypothesis (y-axis)
    pval <- pval(data)

    # Number of reads for the reference
    rReads <- data$seqs$rRead
    
    # Number of reads for the variant
    vReads <- data$seqs$vRead
    
    .plotLODR(data.frame(abund=rReads+vReads, pval=pval, ratio=data$seqs$eAFreq), multiTest=FALSE)
}

plotLODR.LadQuin <- function(data, title=NULL)
{
    # Probability under the null hypothesis (y-axis)
    pval <- pval(data)
    
    # Number of reads for the reference
    abund <- data$seqs$abund

    plotLODR.Plot(data.frame(abund=abund, pval=pval, ratio=data$seqs$expected),
                  title='Significance vs Expression',
                  xname='Expression (log10 Average Reads)',
                  yname='Significance (log10 T-test Stats)',
                  multiTest=FALSE#,
                  #yMin=-2,
                  #yMax=4,
                  #xMin=-3,
                  #xMax=3
                  )
}

plotLODR.TransQuin <- function(data,
                     lvl,
                     title = NULL,
                     chosenFDR = 0.1,
                     xBreaks = NULL,
                     yBreaks = NULL,
                     xLabels = NULL,
                     band = 'pred',
                     xname = 'Average Counts',
                     yname = 'DE Test P-values',
                     shouldTable = FALSE,
                     shouldBand  = FALSE)
{
    .plotLODR(data.frame(abund=data$seqs$mean, pval=pval(data), ratio=abs(data$seqs$expect)), multiTest=FALSE)
}
