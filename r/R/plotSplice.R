#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

#
# Returns the expected concentration for a sequin. The following levels are supported:
#
#   TransQuin: 'exon', 'isoform' and 'gene'
#

expectAbund <- function(data, ids, lvl, mix='A')
{
    stopifnot(lvl == 'exon'    |
                  lvl == 'gene'    |
                  lvl == 'isoform')
    
    stopifnot(class(data) == 'TransQuin' |
                  class(data) == 'VarQuin'   |
                  class(data) == 'MetaQuin'  |
                  class(data) == 'Mixture')
    
    data <- data$mix
    
    switch(lvl, 'gene'    = { data <- data$genes[row.names(data$genes) %in% ids,]       },
           'isoform' = { data <- data$isoforms[row.names(data$isoforms) %in% ids,] },
           'exon'    = { data <- data$exons[row.names(data$exons) %in% ids,]       })
    
    if (is.null(data[mix]))
    {
        error(paste('Unknown mixture:', mix))
    }
    
    return (signif(data[mix][[1]], digits=2))
}

#
# Pool plot is a scatter plot where the variable in the x-axis is categorical. It is a common visualization
# tool for exporing the relationship between sequin groups.
#

expectSplice <- function(data, ids, mix='Mix.A')
{
    stopifnot(class(data) == 'TransQuin')

    splice <- data.frame('prop'=rep(NA, length(ids)))
    row.names(splice) <- ids

    for (gID in row.names(splice))
    {
        # The isoforms for the gene
        isos <- data$mix$isoforms[data$mix$isoforms$GeneID == gID,]
        
        stopifnot(nrow(isos) >= 1)
        
        minor <- isos[which.min(isos$Mix.A),]$Mix.A
        major <- isos[which.max(isos$Mix.A),]$Mix.A

        splice[row.names(splice) == gID,] <- (minor / major)
    }
    
    return (splice)
}


plotSplice <- function(data,
                       xname = 'Log2 expected minor/major',
                       yname = 'Log2 aveagre counts',
                       shouldError = FALSE,
                       cname='Ratio',
                       mix=loadMixture())
{
    require(ggplot2)
    require(RColorBrewer)
    
    stopifnot(class(data) == 'TransQuin')

    # Count table for the replicates
    #samples <- data.frame(A1=data$seqs$A1, A2=data$seqs$A2, A2=data$seqs$A3)
    
    samples <- data.frame(B1=data$seqs$B1, B2=data$seqs$B2, B3=data$seqs$B3) # TODO: FIX THIS!!!!
    row.names(samples) <- row.names(data$seqs)

    # Maximum for the samples
    sMax <- apply(samples, 1, max)
    
    # Minimium for all the samples
    sMin <- apply(samples, 1, min)

    # Arithmetic average of all the samples
    samples$avg <- rowMeans(samples)

    # Mapping from isoforms to genes
    samples$gID <- isoformsToGenes(data, row.names(samples))

    # We can't plot anything that is undefined
    samples <- samples[!is.na(samples$gID),]

    measured <- data.frame(prop=rep(NA, length(unique(samples$gID))))
    row.names(measured) <- unique(samples$gID)
    
    # Expected splicing (minor/major)
    expected <- expectSplice(data, row.names(measured))
    
    for (gID in row.names(measured))
    {
        # The isoforms for the gene
        isos <- samples[samples$gID == gID,]
     
        stopifnot(nrow(isos) >= 1)
        
        print (isos)
        
        minor <- isos[which.min(isos$avg),]$avg
        major <- isos[which.max(isos$avg),]$avg
        prop  <- minor / major
        
        measured[row.names(measured) == gID,] <- prop
    }

    # Expected abundance for each sequin
    abund <- expectAbund(data, lvl='gene', ids=row.names(measured))
    
    data <- data.frame(eProp=expected$prop, prop=measured$prop, abund=abund)
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
