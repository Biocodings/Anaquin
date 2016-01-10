#
#  Copyright (C) 2015 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

#
# Plot the "Limit of Detection ratio" plot. Refer to the ERCC paper for more details.
#

plotLODR <- function(data,
                     choseFDR = 0.1,
                     xname='Average Counts',
                     yname='DE Test P-values',
                     plotTable=TRUE)
{
    require(grid)
    require(qvalue)
    require(locfit)
    require(ggplot2)
    require(gridExtra)
    
    data <- data$seqs

    #
    # Estimate the Q-value for false discovery rate control.
    #
    #    https://github.com/Bioconductor-mirror/erccdashboard/blob/631b691a51db54cb165ff2de1f9a5e06608347bd/R/geneExprTest.R
    #
    
    data$qval <- qvalue(data$pval)$qvalues
    cutoff    <- max(data$pval[data$qval < choseFDR])

    print(paste('FDR threshold:', cutoff))

    xrange <- c(1, 30245)
    
    # Combine the groups solely based on their magnitudes
    data$ratio = abs(data$ratio)
    #data$ratio = as.factor(abs(data$ratio))

    data$x <- data$counts
    data$y <- data$pval
        
    LODR <- function(pval, mn, cutoff, prob)
    {
        # Eg: 0.02590173 to -1.586671
        cutoff <- log10(cutoff)
        
        # Fit a local regression on the log10 scale
        fit <- locfit(log10(pval)~lp(log10(mn)), maxk=300)
        
        X <- preplot(fit,band="pred",newdata=log10(mn))
        
        # See plot of fit
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
    lodr.resPlot<-NULL; set.seed(1)
    lodr.resLess<-NULL; set.seed(1)
    
    #
    # For each group, fit a local regression. We'll also estimate the confidence interval for each
    # model. Refer to the ERCC paper for more details.
    #
    
    prob <- 0.90
    
    for (i in unique(data$ratio))
    {
        t <- data[data$ratio == i,]

        #
        # Not everything can be fitted by local regression. For instance, too few data points would not
        # be appropriate for the model.
        #
        
        tryCatch (
        {
            fit <- locfit(log10(t$y)~lp(log10(t$x)), maxk=300)
            
            x.new <- seq(min(log10(t$x)), max(log10(t$x)), length.out=100)
            X <- preplot(fit,band="pred",newdata=x.new)
            
            x.new <- 10^x.new
            fitLine = 10^(X$fit)
            fitUpper = 10^(X$fit+qnorm(prob)*X$se.fit)
            fitLower = 10^(X$fit-qnorm(prob)*X$se.fit)
            fitData = data.frame(x.new, fitLine, fitUpper, fitLower, ratio = i)
            lineDat = rbind(lineDat, fitData)
        }, error = function(e)
        {
            warning(paste('Failed to fit a local regression for: ', i))
        })
        
        if (i != 0)
        {
            print(paste('Estmating LODR for LFC', i))

            tryCatch (
            {
                t.res <- LODR(t$y, t$x, cutoff=cutoff, prob=prob)
                t.res[-1]<-signif(10^t.res[-1],2)
                    
                if (t.res[1]>.01)
                {
                    t.res[2]<-Inf
                    t.res[3:4]<-NA
                }
                    
                t.resLess <- t.res
                t.resLess[-1][t.resLess[-1] == signif(min(t$x),2)]<-paste("<", signif(min(t$x),2), sep="")
                t.res[-1][t.res[-1]==signif(min(t$x),2)] <- Inf
                    
                lodr.resPlot <- rbind(lodr.resPlot, c(round(abs(as.numeric(i)), 3), t.res))
                lodr.resLess <- rbind(lodr.resLess, c(round(abs(as.numeric(i)), 3), t.resLess))
            }, error = function(e)
            {
                warning(paste('Failed to estimate LODR for: ', i))                
            })
        }
    }
    
    colnames(lodr.resLess)[1:3] <- c("|log2(Fold)|","MinError","Estimate")
    colnames(lodr.resPlot)[1:3] <- c("ratio","MinError","Estimate")
    colnames(lodr.resLess)[1:3] <- c("ratio","MinError","Estimate")
    lodr.resPlot <- as.data.frame(lodr.resPlot)
    lodr.resLess <- as.data.frame(lodr.resLess)

    arrowDat <- data.frame(ratio = lodr.resPlot$ratio, x = lodr.resPlot[,3], y = cutoff, xend = lodr.resPlot[, 3], yend = 0)
    arrowDat$x[grep("<", lodr.resLess[, 3])] <- Inf
    arrowDat <- arrowDat[which(is.finite(arrowDat$x)), ]
    arrowDat$ratio <- as.factor(arrowDat$ratio)
    
    print('Estimation completed')
    
    if (plotTable)
    {
        legendLabels <- c('4', '3', '2', '1')

        lodr.resPlot$ratio <- as.character(legendLabels)
        lodr.resLess$ratio <- as.character(legendLabels) 
        annoTable <- lodr.resLess[-c(2,4,5)]
        
        colnames(annoTable) <- c("ratio", expression("LODR Estimate"))
        cat("\n")
        print(annoTable, quote = FALSE, row.names = FALSE)    
        
        my_table <- tableGrob(d=annoTable, rows=NULL)
    }

    cols <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00") # Colors for the genes
    #cols <- c('#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3', '#fdb462', '#b3de69', '#fccde5', '#d9d9d9', '#bc80bd')

    data$ratio     <- as.factor(data$ratio)
    lineDat$ratio  <- as.factor(lineDat$ratio)
    arrowDat$ratio <- as.factor(arrowDat$ratio)
    
    LODRplot <- ggplot(data, aes(x=x, y=y, colour=ratio)) + 
                             geom_point(size = 6) +
                             xlab(xname) +
                             ylab(yname) +
        
                             scale_x_log10(limits = xrange) + 
                             scale_y_log10(breaks = c(1e-300,1e-200,1e-100,1e-10, 1.00)) +
        
                             geom_ribbon(data  = lineDat, aes(x = x.new, y = fitLine, ymin=fitLower, ymax=fitUpper, fill = ratio),
                                                                alpha = 0.3, colour=NA, show_guide = FALSE) + 
                             geom_line(data = lineDat, aes(x = x.new, y=fitLine, 
                                                                colour = ratio), show_guide = FALSE) +

                             scale_color_manual(values = cols) +
                             scale_fill_manual (values = rev(cols)) +

                             geom_segment(data = arrowDat, 
                                 aes(x = x, y = y, xend = xend, yend = yend, colour = ratio), 
                                     lineend = "round", arrow = grid::arrow(length = grid::unit(0.5, 
                                        "cm")), size = 2, alpha = 0.6) +

                             geom_hline(yintercept = cutoff, linetype = 2, size = 2 ) + # Draw the line for probability threshold
                             theme_bw()

    if (plotTable)
    {
        annotLODRplot <- grid.arrange(arrangeGrob(grobs = list(LODRplot, my_table), ncol = 1, heights = c(2,0.5)))
    }
    else
    {
        print(LODRplot)
    }
    
    lodr <- lodr.resLess[-c(2,4,5)]
    
    #
    # Classify the sequins whether they're above or below the LODR limit.
    #
    
    data$LODR <- NA
    
    for (i in 1:nrow(data))
    {
        if (data$ratio[i] != 0)
        {
            # What's the limit for this sequin?
            limit <- as.numeric(as.character(lodr[lodr$ratio==data$ratio[i],]$Estimate))

            # Classify the sequin based on its average counts (x-axis on the plot)
            data[i,]$LODR <- ifelse(data[i,]$x < limit, 'below', 'above')
        }
    }

    r <- list(lodr.resLess[-c(2,4,5)], data)
    r
}

