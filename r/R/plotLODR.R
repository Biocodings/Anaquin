#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

#
# Plot the "Limit of Detection ratio" plot. Refer to the ERCC paper for more details.
#

plotLODR <- function(data,
                     lvl,
                     title = NULL,
                     choseFDR = 0.1,
                     xBreaks = NULL,
                     yBreaks = NULL,
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
    
    stopifnot(class(data) == 'TransQuin')
    stopifnot(lvl == 'gene' | lvl == 'isoform' | lvl == 'exon')

    # Names of all the features
    names <- names(data)
    
    # Average of normalized count values (x-axis)
    baseMean <- baseMean(data)

    # Probability under the null hypothesis (y-axis)
    pval <- pval(data)
        
    # Expected logFold ratio
    eLogLF <- expectedLF(data, lvl=lvl, ids=sequins(data))

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
        
        X <- preplot(fit, band='pred', newdata=log10(mn))
        
        #plot(fit,band="pred",get.data=TRUE,xlim=range(log10(mn)))
        
        find.mn <- function(mn, fit, cutoff, prob)
        {
            X <- preplot(fit, newdata=mn, band="pred")
            (X$fit+qnorm(prob)*X$se.fit-cutoff)^2
        }
        
        rng.mn <- range(log10(mn))
        
        # Search in sections to get first crossing
        segmented.search <- function(fit)
        {
            X <- preplot(fit, newdata=min(log10(mn)), band="pred")
            
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
                    t.lodr<-optimize(find.mn,c(sum(rng.mn*c(1-ppp+.2,ppp-.2)),
                                               sum(rng.mn*c(1-ppp,ppp))),
                                     fit=fit,cutoff=cutoff,prob=prob)
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
            plot(fit, band="pred", get.data=TRUE, main=paste('Local regression for LFC:', i))
            
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

                      scale_x_log10(limits=c(min(data$baseMean), max(data$baseMean)), breaks=c(arrowDat$x, round(max(data$baseMean)))) +

                      # Draw the fitted lines
                      geom_line(data=lineDat, aes(x=x.new, y=fitLine, colour=ratio), show_guide=FALSE) +

                      geom_segment(data=arrowDat, aes(x = x, y = y, xend = xend, yend = yend, colour = ratio), 
                                     lineend = "round", arrow = grid::arrow(length = grid::unit(0.5, 'cm')), size = 2, alpha = 0.6) +

                      labs(colour='Ratio') +
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
    
    if (shouldBand)
    {
        p <- p + geom_ribbon(data=lineDat, aes(x=x.new, y=fitLine, ymin=fitLower, ymax=fitUpper, fill=ratio),
                             alpha=0.3, colour=NA, show_guide=FALSE)
    }

    if (shouldTable)
    {
        annotLODRplot <- grid.arrange(arrangeGrob(grobs = list(p, my_table), ncol = 1, heights = c(2,0.5)))
    }
    
    print(p)

    #
    # Classify the sequins whether they're above or below the LODR limit.
    #
    lodr <- lodr.resLess[-c(2,4,5)]

#    data$LODR <- NA
    
#    for (i in 1:nrow(data))
#    {
#        if (shouldTable & data$ratio[i] != 0)
#        {
            # What's the limit for this sequin?
#            limit <- as.numeric(as.character(lodr[lodr$ratio==data$ratio[i],]$Estimate))

            # Classify the sequin based on its average counts (x-axis on the plot)
#            data[i,]$LODR <- ifelse(data[i,]$x < limit, 'below', 'above')
#        }
#    }

    return (list(lodr.resLess[-c(2,4,5)], data))
}
