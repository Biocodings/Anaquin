#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#


mLogFSE <- function(data, ids)
{
    stopifnot(class(data) == 'TransQuin' |
                  class(data) == 'VarQuin'   |
                  class(data) == 'MetaQuin')
    
    # Internal representation
    data <- data$seqs
    
    if (is.null(data$lfcSE))
    {
        # TODO: Implement me
    }
    
    return (data$lfcSE)
}


sequins <- function(data)
{
    stopifnot(class(data) == 'TransQuin' |
                  class(data) == 'VarQuin'   |
                  class(data) == 'MetaQuin')
    
    return (row.names(filter(data, 'seqs')))
}

#
# Returns the measured logFold
#
mLogF <- function(data, ids)
{
    stopifnot(class(data) == 'TransQuin' |
                  class(data) == 'VarQuin'   |
                  class(data) == 'MetaQuin')
    
    if (is.null(data$seqs$lfc))
    {
        # TODO: Implement me
    }
    
    return (data$seqs$lfc)
}

#
# Draw an MA plot for differential analysis. The plot depicts the relationship between abundance and variability.
#
#   - x-axis is the measured abundance
#   - y-axis is the measured log-fold change
#
# The following information is required:
#
#   - Sequins
#   - Measured abundance
#   - Measuted log-fold
#   - Expected log-fold
#
plotMA <- function(data,
                   lvl         = 'gene',
                   size        = 4,
                   alpha       = 0.8,
                   title       = NULL,
                   drawHLine   = 'zero',
                   shouldEndo  = TRUE,
                   shouldError = FALSE,
                   shouldSymm  = FALSE,
                   qCutoff     = 0.1,
                   lodr        = NULL,
                   minX        = NULL,
                   maxX        = NULL,
                   minY        = NULL,
                   maxY        = NULL,
                   xname       = 'Log2 Average of Normalized Counts',
                   yname       = 'Log2 Ratio of Normalized Counts')
{
    require(grid)
    require(qvalue)
    require(ggplot2)
    require(gridExtra)

    stopifnot(class(data) == 'TransQuin' |
              class(data) == 'LadQuin')
    stopifnot(lvl == 'gene' | lvl == 'isoform' | lvl == 'exon')

    # Names of all the features
    names <- seqs(data)

    # Abundance (eg: avverage of normalized count values) for the x-axis
    baseMean <- log2(data$seqs$abund)

    # Measured logFold ratio (y-axis)
    logLF <- data$seqs$lfc

    # Expected logFold ratio
    eLogLF <- data$seqs$elfc  #expectLF(data, lvl=lvl, ids=sequins(data))

    # Probability under null hypothesis (optional)
    pvals <- data$seqs$pval

    if (shouldError)
    {
        logFSE <- mLogFSE(data)
    }
    
    data <- data.frame(baseMean=baseMean, logLF=logLF, eLogLF=eLogLF)
    row.names(data) <- names

    if (!is.null(lodr))
    {
        data$LODR <- NA

        for (i in row.names(lodr))
        {
            data[row.names(data)==i,]$LODR <- lodr[row.names(lodr)==i,]
        }
    }
    
    if (!is.null(pvals))
    {
        data$pvals <- pvals
    }

    if (shouldError)
    {
        data$logFSE <- logFSE
    }
    
    data <- data[!is.na(data$baseMean),]
    data <- data[!is.infinite(data$baseMean),]

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

    #
    # Calculating the aspect ratio, and try to maintain the same aspect ratio in both dimension.
    #
    
    if (is.null(minX) || is.null(maxX))
    {
        minX <- round(min(data$baseMean)-1)
        maxX <- round(max(data$baseMean)+1)
    }
    
    if (is.null(minY) || is.null(maxY))
    {
        minY <- round(min(data$logLF)-1)
        maxY <- round(max(data$logLF)+1)
    }
    
    print(minX)
    print(maxX)
    print(minY)
    print(maxY)
    
    xLen <- length(c(minX:maxX))

    xrange <- c(minX, maxX)
    yrange <- c(minY, maxY)

    #
    # Construct the MA plot. We'll draw endogenous features before sequins to avoid overlapping. 
    #

    p <- ggplot(seqs, aes(x=baseMean, y=logLF))                                                
    
    if (shouldEndo)
    {
        endos <- data[is.na(data$eLogLF),]
        
        p <- p + geom_point(data=endos, aes(x=baseMean, y=logLF), colour="grey80", alpha=0.4)
        
        if (!is.null(endos$pval))
        {
            # Convert the p-values into q-values
            endos$qval <- qvalue(endos$pval)$qvalues
            
            p <- p + geom_point(data = endos[endos$qval <= qCutoff,], aes(x=baseMean, y=logLF), colour='pink', alpha=0.20)
        }
    }
    
    p <- p + geom_point(aes(colour=eLogLF), size=size, alpha = alpha)  +
                                                          xlab(xname)  +
                                                          ylab(yname)  +
             coord_cartesian(xlim=xrange, ylim=yrange)                 +
             theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
             scale_y_continuous(breaks=seq(-10, 10, 1))                                +
             labs(colour='Log-Fold')                                                   +
             theme(legend.justification=c(1,0), legend.position=c(1,0)) +
             theme_bw()

    if (!is.null(title))
    {
        p <- p + ggtitle(title)
    }
    
    if (!is.null(lineDat))
    {
        p <- p + geom_hline(data=lineDat, aes(yintercept=logLF, colour=ratio), size=0.5, linetype='longdash')
    }

    if (shouldError)
    {
        p <- p + geom_errorbar(aes(ymax=logLF+logFSE, ymin=logLF-logFSE, colour=eLogLF), size=0.5, alpha=alpha)
    }
    
    if (!is.null(lodr))
    {
        p <- p + geom_point(data=subset(data, LODR=='below'), colour='white', size=2.5)
    }
    
    print(p)

    return (list(xname = xname, yname = yname))
}