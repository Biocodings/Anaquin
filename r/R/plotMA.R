#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

#
# Draw an MA plot for differential analysis. The plot depicts the relationship between abundance and variability.
#
#   - x-axis is the average counts for all samples
#   - y-axis is the measured log-fold change
#

plotMA <- function(data,
                   lvl         = 'gene',
                   size        = 4,
                   alpha       = 0.8,
                   qCutoff     = 0.1,
                   drawHLine   = 'zero',
                   shouldEndo  = TRUE,
                   shouldError = FALSE,
                   shouldSymm  = FALSE,
                   xname       = 'Log2 Average of Normalized Counts',
                   yname       = 'Log2 Ratio of Normalized Counts',
                   shouldLODR=FALSE)
{
    require(grid)
    require(ggplot2)
    require(gridExtra)

    stopifnot(class(data) == 'TransQuin')
    stopifnot(lvl == 'gene' | lvl == 'isoform' | lvl == 'exon')
    
    if (shouldLODR)
    {
        data$status <- NA
        
        #
        # Eg:
        #      R1_43  2705.086868  above
        #      R1_52     7.499938  below
        #
        status <- cutoffs[[2]][c(1,4)]
        status <- status[!is.na(status$LODR),]
        
        for (i in 1:nrow(status))
        {
            t <- status[i,]
            
            if (nrow(data[data$Feature == row.names(t),]) > 0)
            {
                data[data$Feature == row.names(t),]$status <- t$LODR
            }
        }
    }

 #   stats <- function(x, c1, c2)
  #  {
   #     c(mean(log2(x[c2])-log2(x[c1])),  # Ratio of normalized counts (y-axis)
    #      sd(log2(x[c2])-log2(x[c1])),    # Standard deviation
     #     log2(mean(x)))                  # Average normalized counts  (x-axis)
    #} 

    # Names of all the features
    names <- names(data)
    
    # Average of normalized count values (x-axis)
    baseMean <- log2(baseMean(data))

    # Measured logFold ratio (y-axis)
    logLF <- mLogF(data)

    # Expected logFold ratio
    eLogLF <- expectedLF(data, lvl=lvl, ids=sequins(data))

    if (shouldError)
    {
        logFSE <- log2(mLogFSE(data))
    }
    
    data <- data.frame(baseMean=baseMean, logLF=logLF, eLogLF=eLogLF)
    row.names(data) <- names

    if (shouldError)
    {
        data$logFSE <- logFSE
    }
    
    data <- data[!is.na(data$baseMean),]
    data <- data[!is.infinite(data$baseMean),]

    #
    # Eg: { 2.99999, 3.999999999, 1.999999999 }. This can happen if the user supplies like values calculated
    #     in Excel.
    #
    
    if (!all.equal(data$eLogLF, as.integer(data$eLogLF)))
    {
        warning('Non-integer expected log-fold detected. They are now rounded.')
        data$eLogLF <- round(data$eLogLF)
    }
    
    if (shouldSymm)
    {
        data$eLogLF <- abs(data$eLogLF)
    }
    
    # There're only finite possibilities    
    data$eLogLF <- as.factor(data$eLogLF)

    seqs <- data[!is.na(data$eLogLF),]

    #
    # What horizontal lines do we want to draw? Should we draw a line for all groups?
    #
    
    lineDat <- NULL
    
    if (drawHLine == 'zero')
    {
        lineDat <- data.frame(logLF=c(0), ratio=as.factor(c(0)))
    }

    # Number of samples    
#    totCol <- ncol(data)
    
 #   statsDat <- data.frame(t(apply(data[c(1:6)], 1, stats,
  #                                  c1 = c(1:(totCol/2)),
   #                                 c2 = c(((totCol/2)+1):totCol))))
#    colnames(statsDat) <- c("M.Ave", 'M.SD', 'A')

#    data <- cbind(data, statsDat)
#    data <- data[which(is.finite(data$M.Ave)),]

    #
    # Calculating the aspect ratio, and try to maintain the same aspect ratio in both dimension.
    #
    
    minX <- round(min(data$baseMean)-1)
    maxX <- round(max(data$baseMean)+1)
    minY <- round(min(data$logLF)-1)
    maxY <- round(max(data$logLF)+1)
    
    xLen <- length(c(minX:maxX))

    xrange <- c(minX, maxX)
    #yrange <- c(-(xLen/2), (xLen/2))
    yrange <- c(minY, maxY)

    #
    # Construct the MA plot. We'll draw endogenous features before sequins to avoid overlapping. 
    #

    p <- ggplot(seqs, aes(x=baseMean, y=logLF))                                                
    
    if (shouldEndo)
    {
        endos <- data[is.na(data$eLogLF),]
        
        p <- p + geom_point(data=endos, aes(x=baseMean, y=logLF), colour="grey80", alpha=0.5)# +
        #                 geom_point(data = endo[endo$A <= 5,], aes(x = A, y = M.Ave), colour='pink', alpha=0.5)
    }
    
    p <- p + geom_point(aes(colour=eLogLF), size=size, alpha = alpha)  +
                                                          xlab(xname)  +
                                                          ylab(yname)  +
             coord_cartesian(xlim=xrange, ylim=yrange)                 +
             theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
             scale_y_continuous(breaks=seq(-10, 10, 1))                                +
             labs(colour='Log-Ratio')                                       +
             theme(legend.justification=c(1,0), legend.position=c(1,0)) +
             theme_bw()

    if (!is.null(lineDat))
    {
        p <- p + geom_hline(data=lineDat, aes(yintercept=logLF, colour=ratio), size=0.5, linetype='longdash')
    }

    if (shouldError)
    {
        p <- p + geom_errorbar(aes(ymax=logLF+logFSE, ymin=logLF-logFSE, colour=eLogLF), size=0.5, alpha=alpha)
    }
    
    if (shouldLODR)
    {
        p <- p + geom_point(data = subset(data, status == 'below'), colour='white', size=2.5)
    }
    
    print(p)

    return (list(xname = xname, yname = yname))
}