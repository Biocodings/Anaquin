#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

#
# Pool plot is a scatter plot where the variable in the x-axis is categorical. It is a common visualization
# tool for exporing the relationship between sequin groups.
#

plotSplice <- function(data, lvl,
                       xname = '',
                       yname = '',
                       shouldError = FALSE,
                       cname='Ratio',
                       mix=loadMixture())
{
    require(ggplot2)
    require(RColorBrewer)
    
    stopifnot(class(data) == 'TransQuin')
    stopifnot(lvl == 'gene' | lvl == 'isoform' | lvl == 'exon')

    # Count table for the replicates
    samples <- data.frame(A1=data$seqs$A1, A2=data$seqs$A2, A2=data$seqs$A3)
    row.names(samples) <- row.names(data$seqs)
    
    #
    # 1. Calculating the expectation, minimum and maximum
    #

    # Arithemtic average for the samples
    sAve <- rowMeans(samples)

    # Maximum for the samples
    sMax <- apply(samples, 1, max)

    # Minimium for all the samples
    sMin <- apply(samples, 1, min)

    # Expected proportion (eg: minor/major)
    eProp <- data$seqs$prop
    
    # Expected abundance for each sequin
    abund <- expectAbund(data, lvl=lvl, ids=row.names(samples))
    
    data <- data.frame(eProp=eProp, sAve=sAve, sMax=sMax, sMin=sMin, abund=abund)
    
    #seqs <- seqs[seqs$ave!=0,]
    
    # Since the expected abundance varies quite a lot, it's be easier to work on the logarithm scale
    data$abund <- log2(data$abund)
    
    data <- data[data$sAve != 0,]
    
    pal  <- colorRampPalette(brewer.pal(11, "Spectral"))
    cols <- pal(length(unique(data$abund)))
    
    y_min <- round(min(log2(data$sAve)-0.5))
    y_max <- round(max(log2(data$sAve)+0.5))
    
    #    p <- ggplot(data = seqs, aes(x=log2(X), y=ave, colour=as.factor(abund))) +
    p <- ggplot(data=data, aes(x=log2(eProp), y=log2(sAve), colour=abund)) +
                                                        geom_point(size=3) +
                                                               xlab(xname) +
                                                               ylab(yname) +
                                                        labs(colour=cname) +

        #    scale_color_manual(values=cols)      +
        #   scale_fill_manual (values=rev(cols)) +
        scale_colour_gradientn(colours=cols, limits=c(min(data$abund), max(data$abund))) +
        
        scale_x_continuous(breaks=-5:0, labels=c('1/32','1/16','1/8','1/4','1/2','1'), limits=c(-5,0)) +
        scale_y_continuous(breaks=-5:0, labels=c('1/32','1/16','1/8','1/4','1/2','1'), limits=c(-5,0)) +
        
        theme_bw()
    
    if (shouldError)
    {
        p <- p + geom_errorbar(aes(ymax=log2(max), ymin=log2(min)), size=0.3, alpha=0.7)
    }
    
    print(p)
    
    return (list('xname' = xname, 'yname' = yname))
}





