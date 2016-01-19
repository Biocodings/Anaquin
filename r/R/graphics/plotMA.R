#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

#
# Draw an MA plot for differential analysis. The x-axis would be the average counts for all replicates.
# The y-axis would be the measured log-fold change.
#

plotMA <- function(anaquin,
                   lvl='gene',
                   alpha   = 0.8,
                   qCutoff = 0.1,
                   shouldEndo=FALSE,
                   shouldError=FALSE,
                   xname = 'Log2 Average of Normalized Counts',
                   yname = 'Log2 Ratio of Normalized Counts',
                   shouldLODR=FALSE)
{
    require(grid)
    require(ggplot2)
    require(gridExtra)

    stopifnot(class(anaquin) == 'TransQuin')

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

    stats <- function(x, c1, c2)
    {
        c(mean(log2(x[c2])-log2(x[c1])),  # Ratio of normalized counts (y-axis)
          sd(log2(x[c2])-log2(x[c1])),    # Standard deviation
          log2(mean(x)))                  # Average normalized counts  (x-axis)
    } 

    #
    #                     A1  A2  A3  B1  B2  B3
    # ENSG00000000003.14   0   0   0   0   0   0
    # ENSG00000000005.5    0   0   1   0   0   0
    # ENSG00000000419.12 150 115 221  71  77  98
    #
    data = anaquin$seqs
    
    totCol <- ncol(data)

    #
    # Calculating the average counts for each condition
    #

    statsDat <- data.frame(t(apply(data[c(1:6)], 1, stats,
                                    c1 = c(1:(totCol/2)),
                                    c2 = c(((totCol/2)+1):totCol))))
    colnames(statsDat) <- c("M.Ave", 'M.SD', 'A')
    
    data <- cbind(data, statsDat)
    data <- data[which(is.finite(data$M.Ave)),]
    
    #
    # Working out the expected change for each feature.
    #
    
         if (lvl == 'gene')    { si <- row.names(data) %in% row.names(anaquin$mix$genes)    }
    else if (lvl == 'isoform') { si <- row.names(data) %in% row.names(anaquin$mix$isoforms) }
    else if (lvl == 'exon')
    { 
        si <- !is.na(data$ratio)
        #si <- row.names(data) %in% row.names(mix$exons)
    }
    
    #
    # If the expected ratios are not provided, we'll need to calculate from the mixture.
    #
    
    if (is.null(data$ratio))
    {
        data$ratio <- NA
        
        by(data[si,], 1:nrow(data[si,]), function(x)
        {
            # What's the expected log-fold?
            logFold <- expectLF(anaquin, row.names(x), lvl=lvl)
            
            # For this plot, it's okay to plot only the magnitude
            data[si,][row.names(data[si,]) == row.names(x),]$ratio <<- abs(logFold)
        })
    }
    
    data$ratio <- as.factor(as.character(data$ratio))

    seqs <- data[si,]
    stopifnot(sum(is.na(seqs$ratio))  == 0)

    # What lines to draw horizontally?
    lineDat <- data.frame(logFC=c(0), ratio=as.factor(c(0)))

    #
    # Calculating the aspect ratio, and always try to maintain it as a square.
    #
    
    minX <- round(min(seqs$A)-0.5)
    maxX <- round(max(seqs$A)+0.5)
    xLen <- length(c(minX:maxX))

    xrange <- c(minX, maxX)
    yrange <- c(-xLen/2, xLen/2)

    p <- ggplot(seqs, aes(x=A, y=M.Ave))                                                
    
    # Make sure our data points will drawn on top of the endo bins
    if (shouldEndo)
    {
        endo <- data[!si,]
        stopifnot(sum(!is.na(endo$ratio)) == 0)
        
        p <- p + geom_point(data = endo, aes(x = A, y = M.Ave), colour = "grey80", alpha=0.5)# +
        #                 geom_point(data = endo[endo$A <= 5,], aes(x = A, y = M.Ave), colour='pink', alpha=0.5)
    }
    
    p <- p + geom_point(aes(colour=ratio), size=2, alpha = alpha)      +
             xlab(xname)                                               +
             ylab(yname)                                               +
             coord_cartesian(xlim=xrange, ylim = c(-10, 10))           +
                      theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
                      scale_y_continuous(breaks=seq(-10, 10, 1))                                +
             geom_hline(data=lineDat, aes(yintercept=logFC, colour=ratio), size=0.5, linetype='longdash') +
             labs(colour='Ratio')                                       +
             theme(legend.justification=c(1,0), legend.position=c(1,0)) +
             theme_bw()

    if (shouldError)
    {
        p <- p + geom_errorbar(aes(ymax = M.Ave + M.SD, ymin = M.Ave - M.SD,
                               colour=ratio), size = 1, alpha = alphaPoint)
    }
    
    if (shouldLODR)
    {
        p <- p + geom_point(data = subset(data, status == 'below'), colour='white', size=2.5)
    }
    
    print(p)

    return (list(xname = xname, yname = yname))
}