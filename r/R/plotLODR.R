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
    
    for (i in 1:nrow(d))
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

#
# Test at the gene level
#
data <- aqdata(seqs   = c('R1_43','R1_52','R1_91','R2_1','R2_26','R2_32','R2_60','R1_102','R1_11','R1_22','R1_24','R1_32','R1_51','R1_72','R2_67','R2_68','R1_31','R1_93','R2_105','R2_143','R2_47','R2_6','R2_76','R1_101','R1_42','R1_61','R1_73','R2_150','R2_19','R2_27','R2_57','R1_63','R1_83','R2_115','R2_28','R2_38','R1_62','R2_37','R2_41','R2_55','R2_66','R2_72','R2_14','R1_12','R1_53','R2_42','R2_46','R2_73','R1_103','R1_14','R1_33','R1_41','R1_92','R2_153','R2_154','R1_13','R1_21','R1_23','R1_71','R1_81','R1_82','R2_116','R2_117','R2_140','R2_151','R2_152','R2_18','R2_20','R2_24','R2_45','R2_53','R2_54','R2_59','R2_63','R2_65','R2_7','R2_71'),
               counts = c(2705.086868,7.49993785,1.592708121,1.230225393,30245.16372,5.815765526,57.31422832,7.861308236,49.19302776,50.41496679,651.9783904,7.985523155,133.7062386,1.915751089,1.073442671,29.01734957,1137.120055,119.4224592,0.18917613,1.097868629,885.0918089,56.41690991,1.682529281,3.924674515,4603.840008,1.063762344,1009.656343,799.0494006,2483.998333,4.952578082,4.375024791,15084.0071,72.92686167,371.1565224,2.953588894,1.68953155,1.01873706,1.777813952,212.6304813,21327.89169,14764.17376,2.492149405,19030.73113,8.859330357,14.00238258,2.767311968,1.93160761,7.175418831,691.7270976,197.5698542,2.636168967,2381.400168,88.00397943,1.751980522,1690.384465,5658.731489,20826.87313,13.07524022,8130.48873,66.37341747,3234.696422,2.420004433,49.36082846,70.83321355,0.715429424,43.70764867,5519.955865,15.24277683,19.99246051,7.024112983,5.195766319,2147.18809,1.256549113,966.0105647,1.362640676,858.8064128,1.881553695),
               pval   = c(1e-300,0.001169934,0.139450731,0.073978738,9.99999999999997e-311,0.092481602,1.13e-20,0.002833409,2.4e-19,6.48e-21,1.04e-204,0.00057502,1.12e-50,0.69161262,0.45297249,2.21e-13,1.99e-246,2.62e-35,0.581316493,0.358522136,2.7e-247,3.87e-18,0.209233963,0.105514837,1e-305,0.441015205,2.38e-209,2.8e-147,1.09e-295,0.025384555,0.02590173,4.41e-226,7.42e-15,3.42e-83,0.978385972,0.667864534,0.930023203,0.580193256,5e-48,3.07e-122,6.09e-197,0.473146607,7.42e-49,0.307018313,0.475080251,0.922872647,0.714491618,0.475080251,2.29e-26,2.35e-09,0.967411086,3.24e-168,0.012707059,0.575509154,7.88e-31,0.69161262,0.427165219,0.954603427,0.954603427,0.94632463,0.010126637,0.682926366,0.564752562,0.226542104,0.572511389,0.667864534,0.241245248,0.922872647,0.486450354,0.466319074,0.842866512,0.166400526,0.358522136,3.74e-39,0.332977339,0.550500262,0.493058443),
               ratio  = c(4,4,4,4,4,4,4,-4,-4,-4,-4,-4,-4,-4,-4,-4,3,3,3,3,3,3,3,-3,-3,-3,-3,-3,-3,-3,-3,2,2,2,2,2,-2,-2,-2,-2,-2,-2,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
plotLODR(data)

#plotLODR(row.names(m), m$baseMean, m$padj, m$expected.log2FoldChange)
#plotLODR(m$overallMean.FPKM., m$qval, m$expectedLFC, cutoff=0.10)

