#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

.plotLODR <- function(data, choseFDR=0.1, multiTest=TRUE)
{
    stopifnot(!is.null(data$ratio) & !is.null(data$pval) & !is.null(data$abund))
    
    data <- data[!is.na(data$pval),]
    data <- data[!is.na(data$abund),]

    #
    # Estimate the q-value for false discovery rate.
    #
    #    https://github.com/Bioconductor-mirror/erccdashboard/blob/631b691a51db54cb165ff2de1f9a5e06608347bd/R/geneExprTest.R
    #
    
    data$qval <- if(multiTest) qvalue(data$pval)$qvalues else data$pval
    cutoff    <- max(data$pval[data$qval < choseFDR])
    
    print(paste('FDR threshold:', cutoff))
    
    LODR <- function(pval, mn, cutoff, prob)
    {
        # Eg: 0.02590173 to -1.586671
        cutoff <- log10(cutoff)
        
        # Fit a local regression on the log10 scale
        fit <- locfit(log10(pval)~lp(log10(mn)), maxk=300)
        
        X <- preplot(fit, band=locBand, newdata=log10(mn))
        
        #plot(fit,band=locBand,get.data=TRUE,xlim=range(log10(mn)))
        
        find.mn <- function(mn, fit, cutoff, prob)
        {
            X <- preplot(fit, newdata=mn, band=locBand)
            (X$fit+qnorm(prob)*X$se.fit-cutoff)^2
        }
        
        rng.mn <- range(log10(mn))
        
        # Search in sections to get first crossing
        segmented.search <- function(fit)
        {
            X <- preplot(fit, newdata=min(log10(mn)), band=locBand)
            
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
            plot(fit, band=locBand, get.data=TRUE, main=paste('Local regression for LFC:', ratio))
            
            #
            # Generate new data points and use those data points for the prediction band
            #
            
            x.new <- seq(min(log10(t$abund)), max(log10(t$abund)), length.out=100)
            X <- preplot(fit, band=locBand, newdata=x.new)
            
            x.new    <- 10^x.new
            fitLine  <- 10^(X$fit)
            fitUpper <- 10^(X$fit + qnorm(prob) * X$se.fit)  # Upper confidence interval
            fitLower <- 10^(X$fit - qnorm(prob) * X$se.fit)  # Lower confidence interval
            
            fitData <- data.frame(x.new, fitLine, fitUpper, fitLower, ratio=ratio)
            lineDat <- rbind(lineDat, fitData)
        }, error = function(e)
        {
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

                lodr.resPlot <- rbind(lodr.resPlot, c(round(abs(as.numeric(i)), 3), t.res))
                lodr.resLess <- rbind(lodr.resLess, c(round(abs(as.numeric(i)), 3), t.resLess))
            }, error = function(e)
            {
                print(paste('Failed to estimate LODR for:', i))                
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
    cols <- colors(length(data$eLogLF))
    
    data$eLogLF    <- as.factor(data$eLogLF)
    lineDat$ratio  <- as.factor(lineDat$ratio)
    arrowDat$ratio <- as.factor(arrowDat$ratio)    
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
    require(ggplot2)
    
    x <- list(...)
    p <- ggplot(data, aes(x=abund, y=pval, colour=ratio)) + geom_point(size=3) + theme_bw()

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
        p <- p + geom_segment(data=x$arrowDat, aes(x = x, y = y, xend = xend, yend = yend, colour = ratio), 
                     lineend = "round", arrow = grid::arrow(length = grid::unit(0.5, 'cm')), size = 2, alpha = 0.6) +
                      labs(colour='Log-Fold') +
        p <- p + geom_hline(yintercept=cutoff, linetype=2, size=2)
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
        p <- p + scale_x_continuous(breaks=x$xBreaks, labels=x$xLabels)
    }
    #else
    #{
        #p <- p + scale_x_log10(limits=c(min(data$abund), max(data$abund)))
        #p <- p + scale_x_log10(limits=c(min(data$abund), max(data$abund)), breaks=c(arrowDat$x, round(max(data$baseMean))))
    #}
    
    if (!is.null(x$yBreaks))
    {
        #p <- p + scale_y_log10(breaks=x$yBreaks)
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
                     choseFDR = 0.1,
                     xBreaks = NULL,
                     yBreaks = NULL,
                     xLabels = NULL,
                     locBand = 'pred',
                     xname = 'Average Counts',
                     yname = 'DE Test P-values',
                     shouldTable = FALSE,
                     shouldBand  = FALSE)
{
    require(grid)
    require(qvalue)
    require(locfit)
    require(ggplot2)
    require(gridExtra)
    
    stopifnot(lvl == 'gene' | lvl == 'isoform' | lvl == 'exon')

    # Names of all the features
    names <- seqs(data)
    
    # Average of normalized count values (x-axis)
    baseMean <- baseMean(data)

    # Probability under the null hypothesis (y-axis)
    pval <- pval(data)
        
    # Expected logFold ratio
    eLogLF <- expectLF(data, lvl=lvl, ids=sequins(data))

    data <- data.frame(baseMean=baseMean, pval=pval, eLogLF=eLogLF)
    row.names(data) <- names
    
    data <- data[!is.na(data$pval),]
    data <- data[!is.na(data$baseMean),]
    data <- data[!is.infinite(data$baseMean),]
    
    # We can only draw sequins on LODR
    data <- data[!is.na(data$eLogLF),]    
    
    # Combine the groups solely based on their magnitudes
    data$eLogLF = abs(data$eLogLF)

    #
    # Estimate the q-value for false discovery rate.
    #
    #    https://github.com/Bioconductor-mirror/erccdashboard/blob/631b691a51db54cb165ff2de1f9a5e06608347bd/R/geneExprTest.R
    #

    data$qval <- qvalue(data$pval)$qvalues
    cutoff    <- max(data$pval[data$qval < choseFDR])

    print(paste('FDR threshold:', cutoff))

    LODR <- function(pval, mn, cutoff, prob)
    {
        # Eg: 0.02590173 to -1.586671
        cutoff <- log10(cutoff)
        
        # Fit a local regression on the log10 scale
        fit <- locfit(log10(pval)~lp(log10(mn)), maxk=300)
        
        X <- preplot(fit, band=locBand, newdata=log10(mn))
        
        #plot(fit,band=locBand,get.data=TRUE,xlim=range(log10(mn)))
        
        find.mn <- function(mn, fit, cutoff, prob)
        {
            X <- preplot(fit, newdata=mn, band=locBand)
            (X$fit+qnorm(prob)*X$se.fit-cutoff)^2
        }
        
        rng.mn <- range(log10(mn))
        
        # Search in sections to get first crossing
        segmented.search <- function(fit)
        {
            X <- preplot(fit, newdata=min(log10(mn)), band=locBand)
            
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
    # For each group, fit a local regression. We'll also estimate the confidence interval for each
    # model. Refer to the ERCC paper for more details.
    #
    
    prob <- 0.90
    
    for (i in unique(data$eLogLF))
    {
        t <- data[data$eLogLF == i,]

        #
        # 1. Fit a local regression, where each point is approximated by a quadratic function bounded
        #    by a fixed width.
        #

        tryCatch (
        {
            print(paste('Estmating LODR for LFC', i))

            #plot(log10(t$baseMean), log10(t$pval))
            
            # Performs a local regression
            fit <- locfit(log10(t$pval)~lp(log10(t$baseMean)), maxk=300)
            
            # Plot how the points are fitted
            #plot(fit, band=locBand, get.data=TRUE, main=paste('Local regression for LFC:', i))
            
            #
            # Generate new data points and use those data points for the prediction band
            #
            
            x.new <- seq(min(log10(t$baseMean)), max(log10(t$baseMean)), length.out=100)
            X <- preplot(fit, band=locBand, newdata=x.new)

            x.new    <- 10^x.new
            fitLine  <- 10^(X$fit)
            fitUpper <- 10^(X$fit + qnorm(prob) * X$se.fit)  # Upper confidence interval
            fitLower <- 10^(X$fit - qnorm(prob) * X$se.fit)  # Lower confidence interval

            fitData <- data.frame(x.new, fitLine, fitUpper, fitLower, ratio=i)
            lineDat <- rbind(lineDat, fitData)
        }, error = function(e)
        {
            print(paste('Failed to fit a local regression for: ', i))
        })
        
        if (i != 0)
        {
            tryCatch (
            {
                t.res <- LODR(t$pval, t$baseMean, cutoff=cutoff, prob=prob)
                t.res[-1]<-signif(10^t.res[-1],2)

                #if (t.res[1]>.01)
                #{
                 #   t.res[2]<-Inf
                  #  t.res[3:4]<-NA
                #}
                    
                t.resLess <- t.res
                t.resLess[-1][t.resLess[-1] == signif(min(t$baseMean),2)] <- paste("<", signif(min(t$baseMean),2), sep="")
                t.res[-1][t.res[-1]==signif(min(t$baseMean),2)] <- Inf
                    
                lodr.resPlot <- rbind(lodr.resPlot, c(round(abs(as.numeric(i)), 3), t.res))
                lodr.resLess <- rbind(lodr.resLess, c(round(abs(as.numeric(i)), 3), t.resLess))
            }, error = function(e)
            {
                print(paste('Failed to estimate LODR for:', i))                
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
        print(annoTable, quote = FALSE, row.names = FALSE)    
        
        my_table <- tableGrob(d=annoTable, rows=NULL)
    }

    # What color scheme should we use?
    cols <- colors(length(data$eLogLF))
    
    data$eLogLF    <- as.factor(data$eLogLF)
    lineDat$ratio  <- as.factor(lineDat$ratio)
    arrowDat$ratio <- as.factor(arrowDat$ratio)

    p <- ggplot(data, aes(x=baseMean, y=pval, colour=eLogLF)) + 
                      geom_point(size=3)                      +
                      xlab(xname)                             +
                      ylab(yname)                             +

                      # Draw the fitted lines
                      geom_line(data=lineDat, aes(x=x.new, y=fitLine, colour=ratio), show_guide=FALSE) +

                      geom_segment(data=arrowDat, aes(x = x, y = y, xend = xend, yend = yend, colour = ratio), 
                                     lineend = "round", arrow = grid::arrow(length = grid::unit(0.5, 'cm')), size = 2, alpha = 0.6) +

                      labs(colour='Log-Fold') +
                      geom_hline(yintercept=cutoff, linetype=2, size=2) + # Line for probability threshold
                      theme_bw()

    if (!is.null(title))
    {
        p <- p + ggtitle(title)
    }
    
    if (!is.null(yBreaks))
    {
        p <- p + scale_y_log10(breaks=yBreaks)
    }
    
    if (!is.null(xBreaks))
    {
        p <- p + scale_x_log10(breaks=xBreaks, labels=xLabels)
    }
    else
    {
        p <- p + scale_x_log10(limits=c(min(data$baseMean), max(data$baseMean)), breaks=c(arrowDat$x, round(max(data$baseMean))))
    }

    if (shouldBand)
    {
        p <- p + geom_ribbon(data=lineDat, aes(x=x.new, y=fitLine, ymin=fitLower, ymax=fitUpper, fill=ratio),
                             alpha=0.3, colour=NA, show_guide=FALSE)
    }

    if (shouldTable)
    {
        annotLODRplot <- grid.arrange(arrangeGrob(grobs=list(p, my_table), ncol=1, heights=c(2,0.5)))
    }
    
    print(p)

    #
    # Classify the sequins whether they're above or below the LODR limit.
    #
    lodr <- lodr.resLess[-c(2,4,5)]

    # We don't need the sequins without LODR estimate anymore
    data <- data[data$eLogLF != 0,]
    
    # LODR classification for each sequin
    data$LODR <- NA
    
    for (eLogLF in unique(data$eLogLF))
    {
        # The LODR estimate
        limit <- as.numeric(as.character(lodr[lodr$ratio==eLogLF,]$Estimate))
        
        if (length(limit))
        {
            # Classify whether it's below or above the LODR
            data[data$eLogLF==eLogLF,]$LODR <- ifelse(data[data$eLogLF==eLogLF,]$baseMean < limit, 'below', 'above')
        }
    }
    
    lodr$Estimate <- as.numeric(as.character(lodr$Estimate))
    lodr <- lodr[order(lodr$Estimate),]
    print(lodr)

    return (list(lodr=lodr, data=data[c(5)]))
}