#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

expectSplice <- function(data, gIDs)
{
    r <- data.frame('prop'=rep(NA, length(gIDs)), row.names=gIDs)

    for (gID in row.names(r))
    {
        # Relevant isoforms
        x <- data[grep(gID, row.names(data)),]
        
        stopifnot(nrow(x) >= 1)
        
        minor <- x[which.min(x$expected),]$expected
        major <- x[which.max(x$expected),]$expected

        r[row.names(r) == gID,] <- (minor / major)
    }
    
    return (r)
}

plotTMajor <- function(data,
                       xname = 'Log2 expected minor/major',
                       yname = 'Log2 aveagre counts',
                       shouldError = FALSE,
                       cname='Ratio')
{
    require(ggplot2)
    require(RColorBrewer)
    
    stopifnot(class(data) == 'TransQuin')

    samples <- data$seqs$measured
    row.names(samples) <- row.names(data$seqs)

    # Maximum for the samples
    sMax <- apply(samples, 1, max)
    
    # Minimium for all the samples
    sMin <- apply(samples, 1, min)

    # Arithmetic average of all the samples
    samples$avg <- rowMeans(samples)

    # Mapping from isoforms to genes
    samples$gID <- isoformsToGenes(row.names(samples))

    # We can't plot anything that is undefined
    samples <- samples[!is.na(samples$gID),]

    measured <- data.frame(prop=rep(NA, length(unique(samples$gID))))
    row.names(measured) <- unique(samples$gID)
    
    # Expected splicing (minor/major)
    expected <- expectSplice(data$seqs, row.names(measured))
    
    for (gID in row.names(measured))
    {
        # The isoforms for the gene
        isos <- samples[samples$gID == gID,]
     
        stopifnot(nrow(isos) >= 1)
        
        minor <- isos[which.min(isos$avg),]$avg
        major <- isos[which.max(isos$avg),]$avg
        prop  <- minor / major
        
        measured[row.names(measured) == gID,] <- prop
    }

    # Expected abundance for each sequin
    abund <- expected  #expectAbund(data, lvl='gene', ids=row.names(measured))
    
    data <- data.frame(eProp=expected$prop, prop=measured$prop, abund=abund$prop)
    row.names(data) <- row.names(measured)
    
    data <- data[data$prop != 0,]
    data <- data[!is.na(data$prop),]

    # Since the expected abundance varies quite a lot, it's be easier to work on the logarithm scale
    data$abund <- log2(data$abund)
    
    pal  <- colorRampPalette(brewer.pal(11, "Spectral"))
    cols <- pal(length(unique(data$abund)))
    
    y_min <- round(min(log2(data$prop)-0.5))
    y_max <- round(max(log2(data$prop)+0.5))
    
    p <- ggplot(data=data, aes(x=log2(eProp), y=log2(prop), colour=abund)) +
                                                        geom_point(size=3) +
                                                               xlab(xname) +
                                                               ylab(yname) +
                                                        labs(colour=cname) +
        #    scale_color_manual(values=cols)     +
        #   scale_fill_manual (values=rev(cols)) +
        scale_colour_gradientn(colours=cols, limits=c(min(data$abund), max(data$abund))) +
        
        scale_x_continuous(breaks=-5:0, labels=c('1/32','1/16','1/8','1/4','1/2','1'), limits=c(-5,0)) +
        scale_y_continuous(breaks=-5:0, labels=c('1/32','1/16','1/8','1/4','1/2','1'), limits=c(-5,0)) +
        theme_bw()
    
    if (shouldError)
    {
        p <- p + geom_errorbar(aes(ymax=log2(y_max), ymin=log2(y_min)), size=0.3, alpha=0.7)
    }
    
    print(p)
    
    return (list('xname' = xname, 'yname' = yname))
}
