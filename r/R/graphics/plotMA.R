#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

plotMA <- function(data,
                   mix=loadMixture(),
                   alphaPoint = 0.8,
                   qCutoff = 0.1,
                   xname = 'Log2 Average of Normalized Counts',
                   yname = 'Log2 Ratio of Normalized Counts',
                   shouldLODR=FALSE)
{
    require(grid)
    require(ggplot2)
    require(gridExtra)
    
    if (shouldLODR)
    {
        data$status <- NA
        
        #
        # Eg:
        #      R1_43  2705.086868 above
        #      R1_52     7.499938 below
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

    maStats <- function(x, c1, c2)
    {
        c(mean(log2(x[c1])-log2(x[c2])),  # Ratio of normalized counts (y-axis)
          sd(log2(x[c1])-log2(x[c2])),    # Standard deviation
          log2(mean(x)))                  # Average normalized counts (x-axis)
    } 

    totCol <- ncol(data)

    maStatDat <- data.frame(t(apply(data[c(1:6)], 1, maStats,
                                    c1 = c(1:(totCol/2)),
                                    c2 = c(((totCol/2)+1):totCol))))
    colnames(maStatDat) <- c("M.Ave", 'M.SD', 'A')
    
    data <- cbind(data, maStatDat, ratio=NA)
    data <- data[which(is.finite(data$M.Ave)),]
    
    data$qval <- rnorm(nrow(data), qCutoff, 0.1) # TODO: Fix this...

    # Index for sequins
    si <- row.names(data) %in% row.names(mix$genes)

    by(data[si,], 1:nrow(data[si,]), function(x)
    {
        logFold <- loadGene(row.names(x), mix)$logFold
        data[si,][row.names(data[si,]) == row.names(x),]$ratio <<- abs(logFold)
    })

    data$ratio <- as.factor(as.character(data$ratio))

    seqs <- data[si,]
    endo <- data[!si,]

    stopifnot(sum(is.na(seqs$ratio)) == 0)
    stopifnot(sum(!is.na(endo$ratio)) == 0)

    lineDat <- data.frame(logFC=c(0), ratio=c(0))
    lineDat$ratio = as.factor(lineDat$ratio)
    
    # So that the legend starts with an upper case...
    seqs$Ratio <- seqs$ratio

    maPlot <- ggplot(seqs, aes(x = A, y = M.Ave))                                              +
                     geom_point(data = endo,
                                aes(x = A, y = M.Ave), colour = "grey80", alpha = 0.5)         +
                     geom_point(data = endo[endo$A <= 5,],
                                aes(x = A, y = M.Ave), colour = "pink", alpha = 0.5)           +
                     geom_point(aes(colour = Ratio), size = 5, alpha = alphaPoint)             +
                     xlab(xname)                                                               +
                     ylab(yname)                                                               +
                     coord_cartesian(xlim = c(-3,17), ylim = c(-10, 10))                       +
                     geom_errorbar(aes(ymax = M.Ave + M.SD, ymin = M.Ave - M.SD,
                                   colour = Ratio), size = 1, alpha = alphaPoint)              +
                     theme(legend.justification = c(1,0), legend.position=c(1,0))              +
                     theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
                     scale_y_continuous(breaks = seq(-10, 10, 1))                              +
                     geom_hline(data=lineDat, aes(yintercept = logFC, colour = ratio), 
                                   size = 1, linetype = "longdash") +        
        
                     theme_bw()
        
    if (shouldLODR)
    {
        maPlot <- maPlot + geom_point(data = subset(data, status == 'below'), colour = 'white', size = 2.5)
    }
    
    print(maPlot)

    return (list(xname = xname, yname = yname))
}